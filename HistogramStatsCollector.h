#ifndef HISTOGRAMSTATSCOLLECTOR_H
#define HISTOGRAMSTATSCOLLECTOR_H

#pragma once

#include "AbstractStatCollector.h"

namespace BamstatsAlive {

	typedef std::map<size_t, unsigned int> _CoverageHistogramT;

	class CoverageHistogramVisitor;

	class HistogramStatsCollector : public AbstractStatCollector {
		protected:
			unsigned int m_mappingQualHist[256];
			std::map<int32_t, unsigned int> m_fragHist;
			std::map<int32_t, unsigned int> m_lengthHist;
			std::map<std::string, unsigned int> m_refAlnHist;
			_CoverageHistogramT m_covHist;
			unsigned int m_covHistLocs;
			unsigned int m_covHistAccumu;

			const unsigned int kCovHistSkipFactor;

			BamTools::PileupEngine * _pileupEngine;
			CoverageHistogramVisitor * _readDepthHistVisitor;

			virtual void processAlignmentImpl(const BamTools::BamAlignment& al, const BamTools::RefVector& refVector);
			virtual void appendJsonImpl(json_t * jsonRootObj);

		public:
			HistogramStatsCollector(unsigned int skipFactor = 0);
			virtual ~HistogramStatsCollector();


	};
}

#endif
