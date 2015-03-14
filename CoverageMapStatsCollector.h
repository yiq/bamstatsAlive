#ifndef COVERAGEMAPSTATSCOLLECTOR_H
#define COVERAGEMAPSTATSCOLLECTOR_H

#pragma once

#include "AbstractStatCollector.h"
#include "GenomicRegionStore.h"

namespace BamstatsAlive {

	class CoverageMapStatsCollector : public AbstractStatCollector {
		public:
			using coverageHistT = std::map<size_t, unsigned int>;

		protected:
			virtual void processAlignmentImpl(const BamTools::BamAlignment& al, const BamTools::RefVector& refVector);
			virtual void appendJsonImpl(json_t * jsonRootObj);

		private:	
			const GenomicRegionStore::GenomicRegionT *_currentRegion;
			unsigned int * _regionalCoverageMap;
			const coverageHistT& _existingCoverageHist;
			coverageHistT _coverageHist;
			size_t _coveredLength;

		public:
			CoverageMapStatsCollector(const GenomicRegionStore::GenomicRegionT * currentRegion, const coverageHistT& existingHistogram);

			virtual ~CoverageMapStatsCollector();

			coverageHistT getEffectiveHistogram(unsigned int& totalPos); 
	};
}

#endif
