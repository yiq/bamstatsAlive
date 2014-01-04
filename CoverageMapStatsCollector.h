#ifndef COVERAGEMAPSTATSCOLLECTOR_H
#define COVERAGEMAPSTATSCOLLECTOR_H

#pragma once

#include "AbstractStatCollector.h"

namespace BamstatsAlive {

	class CoverageMapStatsCollector : public AbstractStatCollector {
		protected:
			unsigned int m_readDepth[256];

			unsigned int regionStart;
			unsigned int regionLength;

			std::unique_ptr<BamTools::PileupEngine> pileupEngine;
			BamTools::PileupVisitor * visitor;

			virtual void processAlignmentImpl(const BamTools::BamAlignment& al, const BamTools::RefVector& refVector);
			virtual void appendJsonImpl(json_t * jsonRootObj);


		public:
			CoverageMapStatsCollector(unsigned int start, unsigned int length);
			virtual ~CoverageMapStatsCollector();
	};
}

#endif
