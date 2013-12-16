#ifndef BASICSTATSCOLLECTOR_H
#define BASICSTATSCOLLECTOR_H

#pragma once

#include "AbstractStatCollector.h"
#include <map>
#include <string>

namespace BamstatsAlive {

	static std::string const kTotalReads = "total_reads";
	static std::string const kPairedEndReads = "paired_end_reads";
	static std::string const kProperPairs = "proper_pairs";
	static std::string const kMappedReads = "mapped_reads";
	static std::string const kBothMatesMapped = "both_mates_mapped";
	static std::string const kForwardStrands = "forward_strands";
	static std::string const kReverseStrands = "reverse_strands";
	static std::string const kFirstMates = "first_mates";
	static std::string const kSecondMates = "second_mates";
	static std::string const kSingletons = "singletons";
	static std::string const kFailedQC = "failed_qc";
	static std::string const kDuplicates = "duplicates";

	typedef std::map<std::string, unsigned int> StatMapT;

	class BasicStatsCollector : public AbstractStatCollector {

		protected:
			StatMapT m_stats;

			virtual void processAlignmentImpl(const BamTools::BamAlignment& al, const BamTools::RefVector& refVector);
			virtual void appendJsonImpl(json_t * jsonRootObj);

		public:
			BasicStatsCollector();
	};
}

#endif
