#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include <vector>
#include <iostream>

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

using namespace std;

void ProcessAlignment(const BamTools::BamAlignment& al);
void printStats(void);

int main(int argc, char* argv[]) {

	string filename;

	if (argc < 2) 
		filename = "-";
	else 
		filename = argv[1];

	BamTools::BamReader reader;

	reader.Open(filename);

	if(!reader.IsOpen()) {
		cout<<"{\"status\":\"error\", \"message\":\"Cannot open the specified file\"}"<<endl;
		exit(1);
	}

	BamTools::BamAlignment alignment;
	while(reader.GetNextAlignment(alignment)) {
		ProcessAlignment(alignment);
		if(m_numReads % 1000 == 0)
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
}

void printStats(void) {

	cout<<"{";

	cout<<"\"total_reads\":"<<m_numReads<<", ";
	cout<<"\"mapped_reads\":"<<m_numMapped<<", ";
	cout<<"\"forward_strands\":"<<m_numForwardStrand<<", ";
	cout<<"\"reverse_strand\":"<<m_numReverseStrand<<", ";
	cout<<"\"failed_qc\":"<<m_numFailedQC<<", ";
	cout<<"\"duplicates\":"<<m_numDuplicates<<", ";
	cout<<"\"pairedEnd_reads:\":"<<m_numPaired<<", ";
	cout<<"\"proper_pairs:\":"<<m_numProperPair;

	cout<<"}";

  
    // if ( m_numPaired != 0 ) {
    //     cout << "'Proper-pairs':    " << m_numProperPair << "\t(" << ((float)m_numProperPair/m_numPaired)*100 << "%)" << endl;
    //     cout << "Both pairs mapped: " << m_numBothMatesMapped << "\t(" << ((float)m_numBothMatesMapped/m_numPaired)*100 << "%)" << endl;
    //     cout << "Read 1:            " << m_numFirstMate << endl;
    //     cout << "Read 2:            " << m_numSecondMate << endl;
    //     cout << "Singletons:        " << m_numSingletons << "\t(" << ((float)m_numSingletons/m_numPaired)*100 << "%)" << endl;
    // }
    cout << endl;
}
