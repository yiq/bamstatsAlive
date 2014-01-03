#include "AbstractStatCollector.h"
#include "BasicStatsCollector.h"
#include "HistogramStatsCollector.h"
#include "CoverageMapStatsCollector.h"

#include <memory>

static unsigned int totalReads;
static unsigned int updateRate;

static unsigned int regionStart;
static unsigned int regionLength;


using namespace std;
using namespace BamstatsAlive;

void printStatsJansson(AbstractStatCollector& rootStatCollector);

int main(int argc, char* argv[]) {
	
	string filename;
	updateRate = 1000;

	/* process the parameters */
	
	/* In order to use the -s/-l for coverage map statistics, 
	 * One need to make sure that the reads passed in all came
	 * from the same chromosome. Otherwise, the result is undefined.
	 */

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


	/* open BAM file */

	BamTools::BamReader reader;
	reader.Open(filename);
	if(!reader.IsOpen()) {
		cout<<"{\"status\":\"error\", \"message\":\"Cannot open the specified file\"}"<<endl;
		exit(1);
	}


	/* Construct the statistics collectors */

	// NOTICE: The following codes utilize the new c++11 unique_ptr data type
	//         to automatically manage the life-time of heap based objects.
	//         See: http://www.cplusplus.com/reference/memory/unique_ptr/

	// Basic Scalar Statistics
	unique_ptr<BasicStatsCollector> bsc(new BasicStatsCollector());

	// Histogram Statistics
	unique_ptr<HistogramStatsCollector> hsc(new HistogramStatsCollector());
	bsc->addChild(hsc.get());

	// Coverage Map Statistics, only when regionLength is greater than 0
	unique_ptr<CoverageMapStatsCollector> csc;
	if(regionLength != 0) {
		csc.reset(new CoverageMapStatsCollector(regionStart, regionLength));
		bsc->addChild(csc.get());
	}

	/* Process read alignments */

	BamTools::BamAlignment alignment;
	const BamTools::RefVector refVector = reader.GetReferenceData();
	while(reader.GetNextAlignment(alignment)) {
		totalReads++;
		bsc->processAlignment(alignment, refVector);
		if(totalReads > 0 && totalReads % updateRate == 0)
			printStatsJansson(*bsc);
	}

	printStatsJansson(*bsc);
}

void printStatsJansson(AbstractStatCollector& rootStatCollector) {

	// Create the root object that contains everything
	json_t * j_root = json_object();

	// Let the root object of the collector tree create Json
	rootStatCollector.appendJson(j_root);

	// Dump the json
	cout<<json_dumps(j_root, JSON_COMPACT | JSON_ENSURE_ASCII | JSON_PRESERVE_ORDER)<<endl;
}
