#include "CoverageMapStatsCollector.h"

namespace BamstatsAlive {
	class ReadDepthPileupVisitor : public BamTools::PileupVisitor {
		protected:
			unsigned int regionStart;
			unsigned int regionLength;

		public:
			unsigned int m_baseCoverage[256];

			ReadDepthPileupVisitor(unsigned int start, unsigned int length) : regionStart(start), regionLength(length) {

			}

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
}

using namespace BamstatsAlive;
using namespace std;

CoverageMapStatsCollector::CoverageMapStatsCollector(unsigned int start, unsigned int length) : 
	regionStart(start), regionLength(length)
{
	pileupEngine = new BamTools::PileupEngine;
	visitor = new ReadDepthPileupVisitor(regionStart, regionLength);
	pileupEngine->AddVisitor(visitor);
}

CoverageMapStatsCollector::~CoverageMapStatsCollector() {
	delete(pileupEngine);
	delete(visitor);
}

void CoverageMapStatsCollector::processAlignmentImpl(const BamTools::BamAlignment& al, const BamTools::RefVector& refVector) {
	// Read depth stats
	int32_t pos = al.Position;
	if(pos < regionStart || pos > regionStart + regionLength)
		return;
	unsigned int index = (float)(pos - regionStart) / (float)regionLength * 256;
	if(index >= 256) index=255; //Bound Safaguard

	m_readDepth[index]++;

	pileupEngine->AddAlignment(al);
	pileupEngine->Flush();
	delete pileupEngine;
	pileupEngine = new BamTools::PileupEngine;
	pileupEngine->AddVisitor(visitor);
}

void CoverageMapStatsCollector::appendJsonImpl(json_t * jsonRootObj) {
	json_t * j_base_coverage = json_object();
	for(size_t i=0; i<256; i++) {
		if (dynamic_cast<ReadDepthPileupVisitor *>(visitor)->m_baseCoverage[i] > 0) {
			stringstream labelSS; labelSS << i;
			json_object_set_new(j_base_coverage, labelSS.str().c_str(), json_integer(dynamic_cast<ReadDepthPileupVisitor*>(visitor)->m_baseCoverage[i]));
		}
	}
	json_object_set_new(jsonRootObj, "base_coverage", j_base_coverage);

	json_t * j_read_depth = json_object();
	for(size_t i=0; i<256; i++) {
		if (m_readDepth[i] > 0) {
			stringstream labelSS; labelSS << i;
			json_object_set_new(j_read_depth, labelSS.str().c_str(), json_integer(m_readDepth[i]));
		}
	}
	json_object_set_new(jsonRootObj, "read_depth", j_read_depth);
}

