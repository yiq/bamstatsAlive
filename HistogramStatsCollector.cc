#include "HistogramStatsCollector.h"
#include <cmath>

namespace BamstatsAlive {


	class CoverageHistogramVisitor : public BamTools::PileupVisitor {

		_CoverageHistogramT& _histogram;
		unsigned int& _locs;

		public:

			CoverageHistogramVisitor(_CoverageHistogramT& covHist, unsigned int& locs) 
				: BamTools::PileupVisitor(), _histogram(covHist), _locs(locs) {

				;
			}

			void Visit(const BamTools::PileupPosition& pileupData) {
				size_t s = pileupData.PileupAlignments.size();
				_CoverageHistogramT::iterator iter = _histogram.find(s);

				if(iter == _histogram.end()) 	_histogram[s] = 1;
				else 							_histogram[s]++;
				_locs++;
			}
	};
}

using namespace BamstatsAlive;
using namespace std;

HistogramStatsCollector::HistogramStatsCollector(unsigned int skipFactor) : 
	kCovHistSkipFactor(skipFactor), 
	m_covHistLocs(0), 
	m_covHistAccumu(0) 
{

	_pileupEngine = new BamTools::PileupEngine;
	_readDepthHistVisitor = new CoverageHistogramVisitor(m_covHist, m_covHistLocs);
	_pileupEngine->AddVisitor(_readDepthHistVisitor);

	memset(m_mappingQualHist, 0, sizeof(unsigned int) * 256);

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
			unsigned int frag = al.InsertSize; //al.MatePosition - al.Position;
			if(m_fragHist.find(frag) != m_fragHist.end())
				m_fragHist[frag]++;
			else
				m_fragHist[frag] = 1;
		}                   
	}

	if(m_covHistAccumu == 0 && _pileupEngine) _pileupEngine->AddAlignment(al);
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
      double rounded = floor( (it->second / static_cast<double>(m_covHistLocs)) * 1000) / 1000;
      std::stringstream ss(stringstream::in | stringstream::out);
      ss << rounded;
		json_object_set_new(j_cov_hist, labelSS.str().c_str(), json_string( ss.str().c_str() ));
	}
	json_object_set_new(jsonRootObj, "coverage_hist:", j_cov_hist);
}
