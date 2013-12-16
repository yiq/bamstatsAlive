#ifndef HISTOGRAMSTATSCOLLECTOR_H
#define HISTOGRAMSTATSCOLLECTOR_H

#pragma once

#include "AbstractStatCollector.h"

#include <map>
#include <string>
#include <jansson.h>

namespace BamstatsAlive {

	class HistogramStatsCollector : public AbstractStatCollector {
		protected:
			unsigned int m_mappingQualHist[256];
			std::map<int32_t, unsigned int> m_fragHist;
			std::map<int32_t, unsigned int> m_lengthHist;
			std::map<std::string, unsigned int> m_refAlnHist;

			virtual void processAlignmentImpl(const BamTools::BamAlignment& al, const BamTools::RefVector& refVector);
			virtual void appendJsonImpl(json_t * jsonRootObj);
	};
}

#endif
