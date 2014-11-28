#include "AbstractStatCollector.h"
#include "BasicStatsCollector.h"
#include "HistogramStatsCollector.h"
#include "CoverageMapStatsCollector.h"

static unsigned int totalReads;
static unsigned int updateRate;
static unsigned int firstUpdateRate;
static std::string regionJson;
static bool hasRegionSpec = false;

static unsigned int regionStart;
static unsigned int regionLength;


using namespace std;
using namespace BamstatsAlive;

void printStatsJansson(AbstractStatCollector& rootStatCollector);

int main(int argc, char* argv[]) {
	
	string filename;
	updateRate = 1000;
	firstUpdateRate = 0;

	/* process the parameters */
	
	/* In order to use the -s/-l for coverage map statistics, 
	 * One need to make sure that the reads passed in all came
	 * from the same chromosome. Otherwise, the result is undefined.
	 */

	int ch;
	while((ch = getopt(argc, argv, "u:f:s:l:r:")) != -1) {
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
			case 'f':
				firstUpdateRate = atoi(optarg);
				break;
			case 'r':
				regionJson = std::string(optarg);
				hasRegionSpec = true;
				break;
		}
	}

	argc -= optind;
	argv += optind;

	if (argc == 0) 
		filename = "-";
	else 
		filename = *argv;


	/* open BAM file */

	BamTools::BamReader reader;
	reader.Open(filename);
	if(!reader.IsOpen()) {
		cout<<"{\"status\":\"error\", \"message\":\"Cannot open the specified file\"}"<<endl;
		exit(1);
	}

	const BamTools::RefVector refVector = reader.GetReferenceData();

	map<int32_t, string> chromIDNameMap;
	for(size_t i=0; i<refVector.size(); i++) {
		chromIDNameMap[reader.GetReferenceID(refVector[i].RefName)] = refVector[i].RefName;
	}

	/* Construct the statistics collectors */

	// NOTICE: The following codes utilize the new c++11 unique_ptr data type
	//         to automatically manage the life-time of heap based objects.
	//         See: http://www.cplusplus.com/reference/memory/unique_ptr/

	// Basic Scalar Statistics
	BasicStatsCollector bsc;
	HistogramStatsCollector* hsc = NULL;
	GenomicRegionStore *regionStore = NULL;

	// Histogram Statistics
	if(hasRegionSpec) {
		LOGS<<"Has Region Spec"<<std::endl;
		try {
			regionStore = new GenomicRegionStore(regionJson);
			LOGS<<regionStore->regions().size()<<" Regions specified"<<endl;
			hsc = new HistogramStatsCollector(chromIDNameMap, 10, regionStore);

		}
		catch(...) {
			cout<<"{\"status\":\"error\", \"message\":\"Cannot parse region json string\"}"<<endl;
			exit(1);
		}
	}
	else {
		LOGS<<"Does not have region spec"<<std::endl;
		hsc = new HistogramStatsCollector(chromIDNameMap, 10);
	}
	bsc.addChild(hsc);

	// Coverage Map Statistics, only when regionLength is greater than 0
	CoverageMapStatsCollector * csc = NULL;
	if(regionLength != 0) {
		csc = new CoverageMapStatsCollector(regionStart, regionLength);
		bsc.addChild(csc);
	}

	/* Process read alignments */

	BamTools::BamAlignment alignment;
	while(reader.GetNextAlignment(alignment)) {
		totalReads++;
		bsc.processAlignment(alignment, refVector);
		if((totalReads > 0 && totalReads % updateRate == 0) ||
		   (firstUpdateRate>0 && totalReads >= firstUpdateRate)) {
			
			printStatsJansson(bsc);

			// disable first update after it has been fired.
			if(firstUpdateRate > 0) firstUpdateRate = 0;
		}
	}

	// count for all regions from which no read came
	hsc->flushAllRegion();

	printStatsJansson(bsc);

	if(hsc) delete hsc;
	if(regionStore) delete regionStore;
	if(csc) delete csc;
}

void printStatsJansson(AbstractStatCollector& rootStatCollector) {

	// Create the root object that contains everything
	json_t * j_root = json_object();

	// Let the root object of the collector tree create Json
	rootStatCollector.appendJson(j_root);

	// Dump the json
	cout<<json_dumps(j_root, JSON_COMPACT | JSON_ENSURE_ASCII | JSON_PRESERVE_ORDER)<<endl;
}
