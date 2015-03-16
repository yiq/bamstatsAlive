#include "HistogramStatsCollector.h"
#include <cmath>

using namespace BamstatsAlive;
using namespace std;

HistogramStatsCollector::HistogramStatsCollector(std::map<int32_t, std::string>& chromIDNameMap, unsigned int skipFactor, GenomicRegionStore* regionStore) : 
	_chromIDNameMap(chromIDNameMap),
	kCovHistSkipFactor(skipFactor), 
	m_covHistAccumu(0),
	_currentRegion(nullptr),
	_regionStore(regionStore),
	_coverageCollector(nullptr),
	m_covHistTotalPos(0)
{
	memset(m_mappingQualHist, 0, sizeof(unsigned int) * 256);
	memset(m_baseQualHist, 0, sizeof(unsigned int) * 51);

	m_lengthHist.clear();
	m_fragHist.clear();
}

HistogramStatsCollector::~HistogramStatsCollector() {
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

void HistogramStatsCollector::updateRegionalStats(const BamTools::BamAlignment& al, const BamTools::RefVector& refVector) {

	if(!_regionStore) return;

	decltype(_currentRegion) _thisReadRegion = nullptr;

	// locate the region the current read is in
	const GenomicRegionStore::GenomicRegionT& regionWithStart = _regionStore->locateRegion(_chromIDNameMap[al.RefID].c_str(), al.Position);

	if(&regionWithStart != &GenomicRegionStore::kRegionNotFound()) {
		_thisReadRegion = &regionWithStart;
	}
	else {
		const GenomicRegionStore::GenomicRegionT& regionWithEnd = _regionStore->locateRegion(_chromIDNameMap[al.RefID].c_str(), al.Position + al.Length);
		if(&regionWithEnd != &GenomicRegionStore::kRegionNotFound()) {
			_thisReadRegion = &regionWithEnd;
		}
	}

	if(_currentRegion != _thisReadRegion) {
		if(_currentRegion != nullptr && m_covHistAccumu == 0) {
			// switched outside of pile up region
			// update histogram and delete coverage collector
			m_covHist = _coverageCollector->getEffectiveHistogram(m_covHistTotalPos);
			delete _coverageCollector;
			_coverageCollector = nullptr;
		}
		//m_covHistAccumu++;
		_currentRegion = _thisReadRegion;
	}

	if(m_covHistAccumu >= kCovHistSkipFactor) m_covHistAccumu = 0;

	// not even in a pileup region
	if(m_covHistAccumu != 0) return;
	if(_currentRegion == nullptr) return;

	updateBaseQualityHistogram(al);

	// feed pileup
	//_currentRegionLength = _currentRegion->endPos - _currentRegion->startPos + 1;

	if(_coverageCollector == nullptr) {
		_coverageCollector = new CoverageMapStatsCollector(_currentRegion, m_covHist);
	}
	_coverageCollector->processAlignment(al, refVector);
}

void HistogramStatsCollector::processAlignmentImpl(const BamTools::BamAlignment& al, const BamTools::RefVector& refVector) {

	updateReferenceHistogram(al, refVector);

	updateMappingQualityHistogram(al);

	updateReadLengthHistogram(al);

	updateFragmentSizeHistogram(al);

	updateRegionalStats(al, refVector);
}

void HistogramStatsCollector::appendJsonImpl(json_t *jsonRootObj) {

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

   // coverage histogram
   if(_coverageCollector != nullptr)
	   _coverageCollector->appendJson(jsonRootObj);
   else {
	   json_t * j_cov_hist = json_object();

	   for(auto it = m_covHist.begin(); it != m_covHist.end(); it++) {
		   stringstream labelSS; labelSS << it->first;
		   json_object_set_new(j_cov_hist, labelSS.str().c_str(), json_real( it->second / static_cast<double>(m_covHistTotalPos)));
	   }
	   json_object_set_new(jsonRootObj, "coverage_hist", j_cov_hist);

   }
}
