#ifndef HISTOGRAMSTATSCOLLECTOR_H
#define HISTOGRAMSTATSCOLLECTOR_H

#pragma once

#include "AbstractStatCollector.h"
#include "GenomicRegionStore.h"

namespace BamstatsAlive {

	typedef std::map<size_t, unsigned int> _CoverageHistogramT;

	class HistogramStatsCollector : public AbstractStatCollector {
		protected:
			unsigned int m_mappingQualHist[256];
			unsigned int m_baseQualHist[51];
			std::map<int32_t, unsigned int> m_fragHist;
			std::map<int32_t, unsigned int> m_lengthHist;
			std::map<std::string, unsigned int> m_refAlnHist;
			_CoverageHistogramT m_covHist;
			unsigned int m_covHistLocs;
			unsigned int m_covHistAccumu;
			const unsigned int kCovHistSkipFactor;

			virtual void processAlignmentImpl(const BamTools::BamAlignment& al, const BamTools::RefVector& refVector);
			virtual void appendJsonImpl(json_t * jsonRootObj);

		private:
			// Functions to deal with intermitted base coverage calculation
			const GenomicRegionStore::GenomicRegionT *_currentRegion;
			GenomicRegionStore *_regionStore;
			unsigned int * _regionalCoverageMap;


			std::map<int32_t, std::string>& _chromIDNameMap;
			void startBaseCoverageRegion();
			void endBaseCoverageRegion();

			void updateReferenceHistogram(const BamTools::BamAlignment& al, const BamTools::RefVector& refVector);
			void updateMappingQualityHistogram(const BamTools::BamAlignment& al);
			void updateReadLengthHistogram(const BamTools::BamAlignment& al);
			void updateFragmentSizeHistogram(const BamTools::BamAlignment& al);
			void updateBaseQualityHistogram(const BamTools::BamAlignment& al);
			void updateRegionalStats(const BamTools::BamAlignment& al);

		public:
			HistogramStatsCollector(
					std::map<int32_t, std::string>& chromIDNameMap,
					unsigned int skipFactor = 0, 
					GenomicRegionStore* regionStore = NULL);
			virtual ~HistogramStatsCollector();
			void flushActiveRegion();
	};
}

#endif
