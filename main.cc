#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include <src/utils/bamtools_pileup_engine.h>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <assert.h>

static unsigned int m_numReads;
static unsigned int m_numPaired;
static unsigned int m_numProperPair;
static unsigned int m_numMapped;
static unsigned int m_numBothMatesMapped;
static unsigned int m_numForwardStrand;
static unsigned int m_numReverseStrand;
static unsigned int m_numFirstMate;
static unsigned int m_numSecondMate;
static unsigned int m_numSingletons;
static unsigned int m_numFailedQC;
static unsigned int m_numDuplicates;

static unsigned int m_mappingQualHist[256];
static std::map<int32_t, unsigned int>m_lengthHist;
static unsigned int m_readDepth[256];

static unsigned int updateRate;
static unsigned int regionStart;
static unsigned int regionLength;

using namespace std;

void ProcessAlignment(const BamTools::BamAlignment& al);
void printStats(void);

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

			m_readDepth[index]+=depth;
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
	while(reader.GetNextAlignment(alignment)) {
		ProcessAlignment(alignment);
		if(m_numReads % updateRate == 0)
			printStats();
	}

	printStats();
}


void ProcessAlignment(const BamTools::BamAlignment& al) {
  
    // increment total alignment counter
    ++m_numReads;
    
    // incrememt counters for pairing-independent flags
    if ( al.IsDuplicate() ) ++m_numDuplicates;
    if ( al.IsFailedQC()  ) ++m_numFailedQC;
    if ( al.IsMapped()    ) ++m_numMapped;
    
    // increment strand counters
    if ( al.IsReverseStrand() ) 
        ++m_numReverseStrand;
    else 
        ++m_numForwardStrand;

    m_mappingQualHist[al.MapQuality]++;

    if(m_lengthHist.find(al.Length) != m_lengthHist.end())
    	m_lengthHist[al.Length]++;
    else
    	m_lengthHist[al.Length] = 1;
    
    // if alignment is paired-end
    if ( al.IsPaired() ) {
      
        // increment PE counter
        ++m_numPaired;
      
        // increment first mate/second mate counters
        if ( al.IsFirstMate()  ) ++m_numFirstMate;
        if ( al.IsSecondMate() ) ++m_numSecondMate;
        
        // if alignment is mapped, check mate status
        if ( al.IsMapped() ) {
            // if mate mapped
            if ( al.IsMateMapped() ) 
                ++m_numBothMatesMapped;
            // else singleton
            else 
                ++m_numSingletons;
        }
        
        // check for explicit proper pair flag
        if ( al.IsProperPair() )
            ++m_numProperPair;
    }

	// Inform the pileup engine
	if(pileupEngine != NULL) {
		pileupEngine->AddAlignment(al);
	}
}

void printStats(void) {

	cout<<"{";

	cout<<"\"total_reads\":"<<m_numReads<<", ";
	cout<<"\"mapped_reads\":"<<m_numMapped<<", ";
	cout<<"\"forward_strands\":"<<m_numForwardStrand<<", ";
	cout<<"\"reverse_strand\":"<<m_numReverseStrand<<", ";
	cout<<"\"failed_qc\":"<<m_numFailedQC<<", ";
	cout<<"\"duplicates\":"<<m_numDuplicates<<", ";
	cout<<"\"pairedEnd_reads\":"<<m_numPaired<<", ";
	cout<<"\"proper_pairs\":"<<m_numProperPair<<", ";

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
		cout<<"\"read_depth\":{ ";
		firstComma = false;
		for(size_t i=0; i<256; i++) {
			if (m_readDepth[i] > 0) {
				if (firstComma) cout<<", ";
				cout<<"\""<<i<<"\":"<<m_readDepth[i];
				firstComma = true;
			}
		}
		cout<<"},";
	}

	// Output read length hisogram array
	cout<<"\"length_hist\":{ ";
	firstComma = false;

	for(map<int32_t, unsigned int>::iterator it = m_lengthHist.begin(); it!=m_lengthHist.end(); it++) {
		if (firstComma) cout<<", ";
		firstComma = true;

		std::cout<<"\""<<it->first<<"\":"<<it->second;
	}

	cout<<"}";

	
	

	// Finalizing
	cout<<"}";
	
    cout << endl;
}
