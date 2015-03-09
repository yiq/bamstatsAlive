#include "HistogramStatsCollector.h"
#include <cmath>

using namespace BamstatsAlive;
using namespace std;

HistogramStatsCollector::HistogramStatsCollector(std::map<int32_t, std::string>& chromIDNameMap, unsigned int skipFactor, GenomicRegionStore* regionStore) : 
	_chromIDNameMap(chromIDNameMap),
	kCovHistSkipFactor(skipFactor), 
	m_covHistLocs(0),
	m_covHistAccumu(0),
	_currentRegion(nullptr),
	_regionStore(regionStore),
	_regionalCoverageMap(nullptr)
{
	memset(m_mappingQualHist, 0, sizeof(unsigned int) * 256);
	memset(m_baseQualHist, 0, sizeof(unsigned int) * 51);

	m_lengthHist.clear();
	m_fragHist.clear();
}

HistogramStatsCollector::~HistogramStatsCollector() {
	if(_regionalCoverageMap) delete [] _regionalCoverageMap;
}

void HistogramStatsCollector::updateReferenceHistogram(const BamTools::BamAlignment& al, const BamTools::RefVector& refVector) {
	// increment ref aln counter
	if ( al.RefID != -1) m_refAlnHist[ refVector[al.RefID].RefName ]++;
}

void HistogramStatsCollector::updateMappingQualityHistogram(const BamTools::BamAlignment& al) {
    m_mappingQualHist[al.MapQuality]++;    
}

void HistogramStatsCollector::updateReadLengthHistogram(const BamTools::BamAlignment& al) {
	if(m_lengthHist.find(al.Length) != m_lengthHist.end())
		m_lengthHist[al.Length]++;
	else
		m_lengthHist[al.Length] = 1;
}

void HistogramStatsCollector::updateFragmentSizeHistogram(const BamTools::BamAlignment& al) {
	if ( al.IsPaired() && al.IsMapped() && al.IsMateMapped()) {
		if( al.RefID == al.MateRefID && al.MatePosition > al.Position )  {
			int32_t frag = al.InsertSize; //al.MatePosition - al.Position;
			if(m_fragHist.find(frag) != m_fragHist.end())
				m_fragHist[frag]++;
			else
				m_fragHist[frag] = 1;
		}
	}
}

void HistogramStatsCollector::updateBaseQualityHistogram(const BamTools::BamAlignment& al) {
	const unsigned char *q = (const unsigned char *)al.Qualities.c_str();
	if (q[0] != 0xff) {
		for (int i = 0; i < al.Qualities.length(); ++i) {
			unsigned int qual = (unsigned int)(q[i]) - 33;
			if(qual >50) qual = 50;
			m_baseQualHist[qual]++;
		}
	}
}

void HistogramStatsCollector::updateRegionalStats(const BamTools::BamAlignment& al) {

	if(m_covHistAccumu != 0) return;
	
	updateBaseQualityHistogram(al);

	if(_currentRegion == nullptr && _regionStore) {
		// try to identify the region of this read
		const GenomicRegionStore::GenomicRegionT& regionWithStart = _regionStore->locateRegion(_chromIDNameMap[al.RefID].c_str(), al.Position);
		if(&regionWithStart != &GenomicRegionStore::kRegionNotFound()) {
			_currentRegion = &regionWithStart;
		}
		else {
			const GenomicRegionStore::GenomicRegionT& regionWithEnd = _regionStore->locateRegion(_chromIDNameMap[al.RefID].c_str(), al.Position + al.Length);
			if(&regionWithEnd != &GenomicRegionStore::kRegionNotFound()) _currentRegion = &regionWithEnd;
		}
	}

	if(_currentRegion == nullptr) {
		// read outside of the specified regions
		return;
	}

	if(! _currentRegion->contains(_chromIDNameMap[al.RefID].c_str(), al.Position) && ! _currentRegion->contains(_chromIDNameMap[al.RefID].c_str(), al.Position + al.Length)) {
		return;
	}

	// feed pileup
	auto regionalLength = _currentRegion->endPos - _currentRegion->startPos + 1;

	if(_regionalCoverageMap == nullptr) {
		_regionalCoverageMap = new unsigned int[regionalLength];
		memset(_regionalCoverageMap, 0, sizeof(unsigned int) * regionalLength);
	}

	auto readMappedStartPos = al.Position - _currentRegion->startPos;
	auto readMappedEndPos = al.Position + al.Length - _currentRegion->startPos;

	if(readMappedStartPos < 0) readMappedStartPos = 0;
	if(readMappedEndPos >= regionalLength) readMappedEndPos = regionalLength - 1;

	for(auto i=readMappedStartPos; i <= readMappedEndPos; i++) _regionalCoverageMap[i]++;
}

void HistogramStatsCollector::processAlignmentImpl(const BamTools::BamAlignment& al, const BamTools::RefVector& refVector) {

	updateReferenceHistogram(al, refVector);

	updateMappingQualityHistogram(al);

	updateReadLengthHistogram(al);

	updateFragmentSizeHistogram(al);

	updateRegionalStats(al);
}

void HistogramStatsCollector::flushActiveRegion() {
	if(_currentRegion == nullptr) return;

	auto regionalLength = _currentRegion->endPos - _currentRegion->startPos + 1;


	for(size_t i=0; i<regionalLength; i++) {
		auto cov = _regionalCoverageMap[i];
		auto covEntryIter = m_covHist.find(cov);
		if(covEntryIter != m_covHist.end())
			m_covHist[cov]++;
		else
			m_covHist[cov] = 1;
	}
	free(_regionalCoverageMap);
	_regionalCoverageMap = nullptr;
	m_covHistLocs += regionalLength;
}

void HistogramStatsCollector::appendJsonImpl(json_t *jsonRootObj) {
	if(kCovHistSkipFactor != 0 && _regionStore) {
		if(m_covHistAccumu >= kCovHistSkipFactor) {
			// getting ready for another round of piling up

			m_covHistAccumu = 0;
		}
		else {
			++m_covHistAccumu;
			if(_regionalCoverageMap) {
				// if just finished a round of piling up
				flushActiveRegion();
			}
		}
	} else {
		// If kCovHistSkipFactor == 0, reset m_covHistAccumu so that next region
		// is also counted towards coverage histogram
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
