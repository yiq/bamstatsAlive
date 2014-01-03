#include "HistogramStatsCollector.h"

using namespace BamstatsAlive;
using namespace std;

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
}
