#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include <src/utils/bamtools_pileup_engine.h>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <assert.h>

#include <sstream>
#include <jansson.h>

static std::string const kTotalReads = "total_reads";
static std::string const kPairedEndReads = "paired_end_reads";
static std::string const kProperPairs = "proper_pairs";
static std::string const kMappedReads = "mapped_reads";
static std::string const kBothMatesMapped = "both_mates_mapped";
static std::string const kForwardStrands = "forward_strands";
static std::string const kReverseStrands = "reverse_strands";
static std::string const kFirstMates = "first_mates";
static std::string const kSecondMates = "second_mates";
static std::string const kSingletons = "singletons";
static std::string const kFailedQC = "failed_qc";
static std::string const kDuplicates = "duplicates";

static unsigned int m_mappingQualHist[256];
static std::map<int32_t, unsigned int>m_fragHist;
static std::map<int32_t, unsigned int>m_lengthHist;
static std::map<std::string, unsigned int>m_refAlnHist;
static unsigned int m_baseCoverage[256];
static unsigned int m_readDepth[256];

static unsigned int updateRate;
static unsigned int regionStart;
static unsigned int regionLength;


typedef std::map<std::string, unsigned int> StatMapT;

static StatMapT m_stats;




using namespace std;

void ProcessAlignment(const BamTools::BamAlignment& al, const BamTools::RefVector& refVector);
void printStats(void);
void printStatsJansson(void);


class ReadDepthPileupVisitor : public BamTools::PileupVisitor
{
	public:
		virtual void Visit(const BamTools::PileupPosition& pileupData) {

			//determining bin
			int32_t pos = pileupData.Position;
			if(pos < regionStart || pos > regionStart + regionLength)
				return;
			unsigned int index = (float)(pos - regionStart) / (float)regionLength * 256;
			if(index >= 256) index=255; //Bound Safaguard

			int32_t depth = pileupData.PileupAlignments.size();
			for(size_t i=0; i<pileupData.PileupAlignments.size(); i++) {
				if(pileupData.PileupAlignments[i].IsCurrentDeletion)
					depth--;
			}

			m_baseCoverage[index]+=depth;
		}
};

static BamTools::PileupEngine * pileupEngine;
static ReadDepthPileupVisitor * visitor;

int main(int argc, char* argv[]) {

	string filename;
	updateRate = 1000;


	int ch;
	while((ch = getopt(argc, argv, "u:s:l:")) != -1) {
		switch(ch) {
			case 'u':
				updateRate = atoi(optarg);
				break;
			case 's':
				regionStart = atoi(optarg);
				break;
			case 'l':
				regionLength = atoi(optarg);
				break;
		}
	}

	argc -= optind;
	argv += optind;

	if (argc == 0) 
		filename = "-";
	else 
		filename = *argv;

	m_stats[kTotalReads] = 0;
	m_stats[kMappedReads] = 0;
	m_stats[kForwardStrands] = 0;
	m_stats[kReverseStrands] = 0;
	m_stats[kFailedQC] = 0;
	m_stats[kDuplicates] = 0;
	m_stats[kPairedEndReads] = 0;
	m_stats[kProperPairs] = 0;
	m_stats[kBothMatesMapped] = 0;
	m_stats[kFirstMates] = 0;
	m_stats[kSecondMates] = 0;
	m_stats[kSingletons] = 0;

	BamTools::BamReader reader;

	reader.Open(filename);

	if(!reader.IsOpen()) {
		cout<<"{\"status\":\"error\", \"message\":\"Cannot open the specified file\"}"<<endl;
		exit(1);
	}

	if(regionLength != 0) {
		visitor = new ReadDepthPileupVisitor();
		pileupEngine = new BamTools::PileupEngine;
		pileupEngine->AddVisitor(visitor);
	}

	BamTools::BamAlignment alignment;
	const BamTools::RefVector refVector = reader.GetReferenceData();
	while(reader.GetNextAlignment(alignment)) {
		ProcessAlignment(alignment, refVector);
		if(m_stats[kTotalReads] > 0 && m_stats[kTotalReads] % updateRate == 0)
			printStats();
	}

	//printStats();
	printStatsJansson();
}


void ProcessAlignment(const BamTools::BamAlignment& al, const BamTools::RefVector& refVector) {
  
    // increment total alignment counter
    ++m_stats[kTotalReads];
    
    // increment ref aln counter
    if ( al.RefID != -1) m_refAlnHist[ refVector[al.RefID].RefName ]++;
    
    // incrememt counters for pairing-independent flags
    if ( al.IsDuplicate() ) ++m_stats[kDuplicates];
    if ( al.IsFailedQC()  ) ++m_stats[kFailedQC];
    if ( al.IsMapped()    ) ++m_stats[kMappedReads];
    
    // increment strand counters
    if ( al.IsReverseStrand() ) 
        ++m_stats[kReverseStrands];
    else 
        ++m_stats[kForwardStrands];

    m_mappingQualHist[al.MapQuality]++;    

    if(m_lengthHist.find(al.Length) != m_lengthHist.end())
    	m_lengthHist[al.Length]++;
    else
    	m_lengthHist[al.Length] = 1;
    
    // if alignment is paired-end
    if ( al.IsPaired() ) {
      
        // increment PE counter
        ++m_stats[kPairedEndReads];
      
        // increment first mate/second mate counters
        if ( al.IsFirstMate()  ) ++m_stats[kFirstMates];
        if ( al.IsSecondMate() ) ++m_stats[kSecondMates];
        
        // if alignment is mapped, check mate status
        if ( al.IsMapped() ) {
            // if mate mapped
            if ( al.IsMateMapped() )  {
                ++m_stats[kBothMatesMapped];
                // record fragment length of first pairs
                if( al.MatePosition > al.Position )  {
                   unsigned int frag = al.MatePosition - al.Position;
                   if(m_fragHist.find(frag) != m_fragHist.end())
                    	m_fragHist[frag]++;
                   else
                    	m_fragHist[frag] = 1;
                }                   
             }
            // else singleton
            else 
                ++m_stats[kSingletons];
        }
        
        // check for explicit proper pair flag
        if ( al.IsProperPair() )
            ++m_stats[kProperPairs];
    }

    // Read depth stats
    if (regionLength > 0)
    {
    	int32_t pos = al.Position;
    	if(pos < regionStart || pos > regionStart + regionLength)
    		return;
    	unsigned int index = (float)(pos - regionStart) / (float)regionLength * 256;
		if(index >= 256) index=255; //Bound Safaguard
   
		m_readDepth[index]++;
    }

    // Inform the pileup engine
	if(pileupEngine != NULL) {
		pileupEngine->AddAlignment(al);
		pileupEngine->Flush();
		delete pileupEngine;
		pileupEngine = new BamTools::PileupEngine;
		pileupEngine->AddVisitor(visitor);
	}

}

void printStats(void) {

	cout<<"{";

	StatMapT::iterator iter;
	for(iter = m_stats.begin(); iter != m_stats.end(); iter++) {
		cout<<"\""<<iter->first<<"\":"<<iter->second<<", ";
	}

	// Output mapping quality array
	cout<<"\"mapq_hist\":{ ";
	bool firstComma = false;
	for(size_t i=0; i<256; i++) {
		if (m_mappingQualHist[i] > 0) {
			if (firstComma) cout<<", ";
			cout<<"\""<<i<<"\":"<<m_mappingQualHist[i];
			firstComma = true;
		}
	}
	cout<<"}, ";

	// Output read depth histogram array
	if(regionLength > 0) {
		cout<<"\"base_coverage\":{ ";
		firstComma = false;
		for(size_t i=0; i<256; i++) {
			if (m_baseCoverage[i] > 0) {
				if (firstComma) cout<<", ";
				cout<<"\""<<i<<"\":"<<m_baseCoverage[i];
				firstComma = true;
			}
		}
		cout<<"}, ";

		cout<<"\"read_depth\":{ ";
		firstComma = false;
		for(size_t i=0; i<256; i++) {
			if (m_readDepth[i] > 0) {
				if (firstComma) cout<<", ";
				cout<<"\""<<i<<"\":"<<m_readDepth[i];
				firstComma = true;
			}
		}
		cout<<"}, ";
	}

	// Output read length histogram array
	cout<<"\"length_hist\":{ ";
	firstComma = false;

	for(map<int32_t, unsigned int>::iterator it = m_lengthHist.begin(); it!=m_lengthHist.end(); it++) {
		if (firstComma) cout<<", ";
		firstComma = true;

		std::cout<<"\""<<it->first<<"\":"<<it->second;
	}
	cout<<"}, ";
	
	// Output ref alignment histogram array
	cout<<"\"refAln_hist\":{ ";
	firstComma = false;

	for(map<std::string, unsigned int>::iterator it = m_refAlnHist.begin(); it!=m_refAlnHist.end(); it++) {
		if (firstComma) cout<<", ";
		firstComma = true;

		std::cout<<"\""<<it->first<<"\":"<<it->second;
	}
	cout<<"}, ";
	
	
	// Output fragment length hisogram array
	cout<<"\"frag_hist\":{ ";
	firstComma = false;

	for(map<int32_t, unsigned int>::iterator it = m_fragHist.begin(); it!=m_fragHist.end(); it++) {
		if (firstComma) cout<<", ";
		firstComma = true;

		std::cout<<"\""<<it->first<<"\":"<<it->second;
	}
	cout<<"}";

	
	

	// Finalizing
	cout<<"}";
	
    cout << endl;
}

void printStatsJansson(void) {

	// Create the root object that contains everything
	json_t * j_root = json_object();

	// Scalar Statistics
	StatMapT::iterator iter;
	for(iter = m_stats.begin(); iter != m_stats.end(); iter++) 
		json_object_set_new(j_root, iter->first.c_str(), json_integer(iter->second));

	// Mapping quality map
	json_t * j_mapq_hist = json_object();
	for(size_t i=0; i<256; i++) {
		if (m_mappingQualHist[i] > 0) {
			stringstream labelSS; labelSS << i;
			json_object_set_new(j_mapq_hist, labelSS.str().c_str(), json_integer(m_mappingQualHist[i]));
		}
	}
	json_object_set_new(j_root, "mapq_hist", j_mapq_hist);

	// Output read depth histogram array
	if(regionLength > 0) {
		json_t * j_base_coverage = json_object();
		for(size_t i=0; i<256; i++) {
			if (m_baseCoverage[i] > 0) {
				stringstream labelSS; labelSS << i;
				json_object_set_new(j_base_coverage, labelSS.str().c_str(), json_integer(m_baseCoverage[i]));
			}
		}
		json_object_set_new(j_root, "base_coverage", j_base_coverage);

		json_t * j_read_depth = json_object();
		for(size_t i=0; i<256; i++) {
			if (m_readDepth[i] > 0) {
				stringstream labelSS; labelSS << i;
				json_object_set_new(j_read_depth, labelSS.str().c_str(), json_integer(m_readDepth[i]));
			}
		}
		json_object_set_new(j_root, "read_depth", j_read_depth);
	}

	// Read length histogram array
	json_t * j_length_hist = json_object();
	for(map<int32_t, unsigned int>::iterator it = m_lengthHist.begin(); it!=m_lengthHist.end(); it++) {
		stringstream labelSS; labelSS << it->first;
		json_object_set_new(j_length_hist, labelSS.str().c_str(), json_integer(it->second));
	}
	json_object_set_new(j_root, "length_hist", j_length_hist);
	
	// Reference alignment histogram array
	json_t * j_refAln_hist = json_object();
	for(map<std::string, unsigned int>::iterator it = m_refAlnHist.begin(); it!=m_refAlnHist.end(); it++) {
		stringstream labelSS; labelSS << it->first;
		json_object_set_new(j_refAln_hist, labelSS.str().c_str(), json_integer(it->second));
	}
	json_object_set_new(j_root, "refAln_hist", j_refAln_hist);
	
	// Fragment length hisogram array
	
	json_t * j_frag_hist = json_object();
	for(map<int32_t, unsigned int>::iterator it = m_fragHist.begin(); it!=m_fragHist.end(); it++) {
		stringstream labelSS; labelSS<<it->first;
		json_object_set_new(j_frag_hist, labelSS.str().c_str(), json_integer(it->second));
	}
	json_object_set_new(j_root, "frag_hist", j_frag_hist);

	cout<<json_dumps(j_root, JSON_COMPACT | JSON_ENSURE_ASCII | JSON_PRESERVE_ORDER)<<endl;
}
