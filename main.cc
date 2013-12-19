#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <assert.h>

#include <jansson.h>

#include "AbstractStatCollector.h"
#include "BasicStatsCollector.h"
#include "HistogramStatsCollector.h"
#include "CoverageMapStatsCollector.h"

static unsigned int totalReads;
static unsigned int updateRate;

static unsigned int regionStart;
static unsigned int regionLength;


using namespace std;
using namespace BamstatsAlive;

void printStatsJansson(AbstractStatCollector * rootStatCollector);

/** statistics collectors **/
BasicStatsCollector * bsc;
HistogramStatsCollector * hsc;
CoverageMapStatsCollector * csc;

int main(int argc, char* argv[]) {
	
	string filename;
	updateRate = 1000;

	/* process the parameters */

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

	bsc = new BasicStatsCollector();

	hsc = new HistogramStatsCollector();
	//bsc->addChild(hsc);

	if(regionLength != 0) {
		csc = new BamstatsAlive::CoverageMapStatsCollector(regionStart, regionLength);
		//bsc->addChild(csc);
	}


	/* Process read alignments */

	BamTools::BamAlignment alignment;
	const BamTools::RefVector refVector = reader.GetReferenceData();
	while(reader.GetNextAlignment(alignment)) {
		totalReads++;
		bsc->processAlignment(alignment, refVector);
		if(totalReads > 0 && totalReads % updateRate == 0)
			printStatsJansson(bsc);
	}

	printStatsJansson(bsc);


	/* clean up */

	delete bsc;
	delete hsc;
	
	if(regionLength != 0)
		delete csc;
}

void printStatsJansson(AbstractStatCollector * rootStatCollector) {

	// Create the root object that contains everything
	json_t * j_root = json_object();

	// Let the root object of the collector tree create Json
	rootStatCollector->appendJson(j_root);

	// Dump the json
	cout<<json_dumps(j_root, JSON_COMPACT | JSON_ENSURE_ASCII | JSON_PRESERVE_ORDER)<<endl;
}
