#include "HistogramStatsCollector.h"
#include <cmath>

namespace BamstatsAlive {

	static std::map<const GenomicRegionStore::GenomicRegionT *, bool> hasRegionBeenConsidered;
	static std::map<const GenomicRegionStore::GenomicRegionT *, bool> isRegionContainingReads;

	class CoverageHistogramVisitor : public BamTools::PileupVisitor {

		_CoverageHistogramT& _histogram;
		unsigned int& _locs;
		bool hasSeenFirst;
		unsigned int lastLoc;

		const GenomicRegionStore::GenomicRegionT * _currentRegion;

		public:

			CoverageHistogramVisitor(_CoverageHistogramT& covHist, unsigned int& locs) 
				: BamTools::PileupVisitor(), 
				_histogram(covHist), 
				_locs(locs), 
				lastLoc(0),
				hasSeenFirst(false) {

					;
				}

			void SetRegion(const GenomicRegionStore::GenomicRegionT * region) {
				_currentRegion = region;
				hasSeenFirst = false;
			}

			void Visit(const BamTools::PileupPosition& pileupData) {

				if(_currentRegion != NULL && hasRegionBeenConsidered.find(_currentRegion) == hasRegionBeenConsidered.end()) {
					int32_t prefixGapSize = pileupData.Position - _currentRegion->startPos;
					LOGS<<"First position in a region, with "<<prefixGapSize<<" prefix gap size"<<std::endl;
					if(prefixGapSize >= 0) {
						_histogram[0] += prefixGapSize * 2;
						_locs += prefixGapSize * 2;
					}
					hasRegionBeenConsidered[_currentRegion] = true;
				}

				size_t s = pileupData.PileupAlignments.size();
				lastLoc = pileupData.Position;

				LOGS<<"PileupVisit: "<<pileupData.Position<<", coverage: "<<s<<std::endl;

				_CoverageHistogramT::iterator iter = _histogram.find(s);

				if(iter == _histogram.end()) 	_histogram[s] = 1;
				else 							_histogram[s]++;
				_locs++;
			}
	};
}

using namespace BamstatsAlive;
using namespace std;

HistogramStatsCollector::HistogramStatsCollector(std::map<int32_t, std::string>& chromIDNameMap, unsigned int skipFactor, GenomicRegionStore* regionStore) : 
	_chromIDNameMap(chromIDNameMap),
	kCovHistSkipFactor(skipFactor), 
	m_covHistLocs(0), 
	m_covHistAccumu(0),
	_currentRegion(NULL),
	_regionStore(regionStore)
{

	_pileupEngine = new BamTools::PileupEngine;
	_readDepthHistVisitor = new CoverageHistogramVisitor(m_covHist, m_covHistLocs);
	_pileupEngine->AddVisitor(_readDepthHistVisitor);


	memset(m_mappingQualHist, 0, sizeof(unsigned int) * 256);
	memset(m_baseQualHist, 0, sizeof(unsigned int) * 51);

	m_lengthHist.clear();
	m_fragHist.clear();
}

HistogramStatsCollector::~HistogramStatsCollector() {
	if(_readDepthHistVisitor) delete _readDepthHistVisitor;
	if(_pileupEngine) delete _pileupEngine;
}

void HistogramStatsCollector::processAlignmentImpl(const BamTools::BamAlignment& al, const BamTools::RefVector& refVector) {

    // increment ref aln counter
    if ( al.RefID != -1) m_refAlnHist[ refVector[al.RefID].RefName ]++;
    
    m_mappingQualHist[al.MapQuality]++;    

    if(m_lengthHist.find(al.Length) != m_lengthHist.end())
    	m_lengthHist[al.Length]++;
    else
    	m_lengthHist[al.Length] = 1;

	// if alignment is paired-end
	if ( al.IsPaired() && al.IsMapped() && al.IsMateMapped()) {
		if( al.RefID == al.MateRefID && al.MatePosition > al.Position )  {
			int32_t frag = al.InsertSize; //al.MatePosition - al.Position;
			if(m_fragHist.find(frag) != m_fragHist.end())
				m_fragHist[frag]++;
			else
				m_fragHist[frag] = 1;
		}                   
	}


	const GenomicRegionStore::GenomicRegionT * readRegion = NULL;
	if(_regionStore != NULL) {
		const GenomicRegionStore::GenomicRegionT& readRegionStart = _regionStore->locateRegion(_chromIDNameMap[al.RefID].c_str(), al.Position);
		if(&readRegionStart != &GenomicRegionStore::kRegionNotFound()) {
			readRegion = &readRegionStart;
		}
		else {
			const GenomicRegionStore::GenomicRegionT& readRegionEnd = _regionStore->locateRegion(_chromIDNameMap[al.RefID].c_str(), al.Position + al.Length);
			if(&readRegionEnd != &GenomicRegionStore::kRegionNotFound()) {
				readRegion = &readRegionEnd;
			}
		}
	}

	if(readRegion != NULL) {
		isRegionContainingReads[readRegion] = true;
	}


	if(m_covHistAccumu == 0 && _pileupEngine) {
		// add base quality histogram
		const unsigned char *q = (const unsigned char *)al.Qualities.c_str();
		if (q[0] != 0xff) {
			for (int i = 0; i < al.Qualities.length(); ++i) {
				unsigned int qual = (unsigned int)(q[i]) - 33;
				if(qual >50) qual = 50;
				m_baseQualHist[qual]++;
			}
		}

		if(_currentRegion == NULL && _regionStore) {
			// try to identify the region of this read
			const GenomicRegionStore::GenomicRegionT& regionWithStart = _regionStore->locateRegion(_chromIDNameMap[al.RefID].c_str(), al.Position);
			if(&regionWithStart != &GenomicRegionStore::kRegionNotFound()) {
				_currentRegion = &regionWithStart;
			}
			else {
				const GenomicRegionStore::GenomicRegionT& regionWithEnd = _regionStore->locateRegion(_chromIDNameMap[al.RefID].c_str(), al.Position + al.Length);
				if(&regionWithEnd != &GenomicRegionStore::kRegionNotFound()) _currentRegion = &regionWithEnd;
			}
			_readDepthHistVisitor->SetRegion(_currentRegion);
		}

		if(_currentRegion)
			LOGS<<">> Looking at read "<<al.Name<<", "<<al.Position<<" - "<<al.Position + al.Length<<std::endl;
		else
			LOGS<<"!! Looking at read "<<al.Name<<", "<<al.Position<<" - "<<al.Position + al.Length<<std::endl;

		// check if active region is in place
		if(_currentRegion) {
			if(_currentRegion->contains(_chromIDNameMap[al.RefID].c_str(), al.Position) || _currentRegion->contains(_chromIDNameMap[al.RefID].c_str(), al.Position + al.Length))
				_pileupEngine->AddAlignment(al);
		}
		else {
			// if somehow there is no active region, yet pileupEngine is active, feed the read anyway
			if(_regionStore == NULL) _pileupEngine->AddAlignment(al);
		}
	}
}

void HistogramStatsCollector::flushAllRegion() {
	if(_regionStore) {
		for(std::vector<GenomicRegionStore::GenomicRegionT>::const_iterator it = _regionStore->regions().begin(); it!=_regionStore->regions().end(); it++) {
			LOGS<<"REGION "<<it->startPos<<" -- "<<it->endPos;
			if(isRegionContainingReads.find(&(*it)) != isRegionContainingReads.end()) {
				LOGS <<" has read, skipping" <<std::endl;
				continue;
			}
			LOGS << " does not have read, counting 0s"<<std::endl;

			m_covHist[0] += (it->endPos - it->startPos);
			m_covHistLocs += (it->endPos - it->startPos);
		}
	}
}

void HistogramStatsCollector::appendJsonImpl(json_t *jsonRootObj) {


	if(kCovHistSkipFactor != 0) {
		if(m_covHistAccumu >= kCovHistSkipFactor) {
			// getting ready for another round of piling up

			m_covHistAccumu = 0;
			_pileupEngine = new BamTools::PileupEngine;
			//_readDepthHistVisitor = new CoverageHistogramVisitor(m_covHist, m_covHistLocs);
			_pileupEngine->AddVisitor(_readDepthHistVisitor);
		}
		else {
			++m_covHistAccumu;
			if(_pileupEngine) {
				// if just finished a round of piling up

				_pileupEngine->Flush();
				delete _pileupEngine;
				_pileupEngine = NULL;

				// Disable current region, so that next time pileup engine is on,
				// current region can be recalculated based on new reads
				if(_currentRegion)
					_currentRegion = NULL;
			}
		}
	} else {
		m_covHistAccumu = 0;
	}

   // Mapping quality map
   json_t * j_mapq_hist = json_object();
   for(size_t i=0; i<256; i++) {
      if (m_mappingQualHist[i] > 0) {
         stringstream labelSS; labelSS << i;
         json_object_set_new(j_mapq_hist, labelSS.str().c_str(), json_integer(m_mappingQualHist[i]));
      }
   }
   json_object_set_new(jsonRootObj, "mapq_hist", j_mapq_hist);
   
   // Base quality map
      json_t * j_baseq_hist = json_object();
      for(size_t i=0; i<=50; i++) {
         if (m_baseQualHist[i] > 0) {
            stringstream labelSS; labelSS << i;
            json_object_set_new(j_baseq_hist, labelSS.str().c_str(), json_integer(m_baseQualHist[i]));
         }
      }
      json_object_set_new(jsonRootObj, "baseq_hist", j_baseq_hist);  
   
   // Fragment length hisogram array
   json_t * j_frag_hist = json_object();
   for(map<int32_t, unsigned int>::iterator it = m_fragHist.begin(); it!=m_fragHist.end(); it++) {
      stringstream labelSS; labelSS<<it->first;
      json_object_set_new(j_frag_hist, labelSS.str().c_str(), json_integer(it->second));
   }
   json_object_set_new(jsonRootObj, "frag_hist", j_frag_hist);
   
   // Read length histogram array
   json_t * j_length_hist = json_object();
   for(map<int32_t, unsigned int>::iterator it = m_lengthHist.begin(); it!=m_lengthHist.end(); it++) {
      stringstream labelSS; labelSS << it->first;
      json_object_set_new(j_length_hist, labelSS.str().c_str(), json_integer(it->second));
   }
   json_object_set_new(jsonRootObj, "length_hist", j_length_hist);
   
   // Reference alignment histogram array
   json_t * j_refAln_hist = json_object();
   for(map<std::string, unsigned int>::iterator it = m_refAlnHist.begin(); it!=m_refAlnHist.end(); it++) {
      stringstream labelSS; labelSS << it->first;
      json_object_set_new(j_refAln_hist, labelSS.str().c_str(), json_integer(it->second));
   }
   json_object_set_new(jsonRootObj, "refAln_hist", j_refAln_hist);

	// Coverage Histogram
	json_t * j_cov_hist = json_object();
	for(_CoverageHistogramT::const_iterator it = m_covHist.begin(); it != m_covHist.end(); it++) {
		stringstream labelSS; labelSS << it->first;
		json_object_set_new(j_cov_hist, labelSS.str().c_str(), json_real( it->second / static_cast<double>(m_covHistLocs)));
	}
	json_object_set_new(jsonRootObj, "coverage_hist", j_cov_hist);
}
