#include "BasicStatsCollector.h"
#include <iostream>

using namespace BamstatsAlive;

BasicStatsCollector::BasicStatsCollector() {
	m_stats[kTotalReads] = 0;
	m_stats[kMappedReads] = 0;
	m_stats[kForwardStrands] = 0;
	m_stats[kReverseStrands] = 0;
	m_stats[kFailedQC] = 0;
	m_stats[kDuplicates] = 0;
	m_stats[kPairedEndReads] = 0;
	m_stats[kProperPairs] = 0;
	m_stats[kBothMatesMapped] = 0;
	m_stats[kFirstMates] = 0;
	m_stats[kSecondMates] = 0;
	m_stats[kSingletons] = 0;

	m_stats.clear();

	StatMapT::iterator iter;
	for(iter = m_stats.begin(); iter != m_stats.end(); iter++) {
		std::cerr<<"Initializing: "<<iter->first<<std::endl;
	}

}

void BasicStatsCollector::processAlignmentImpl(const BamTools::BamAlignment& al, const BamTools::RefVector& refVector) {
	// increment total alignment counter
	++m_stats[kTotalReads];

	// incrememt counters for pairing-independent flags
	if ( al.IsDuplicate() ) ++m_stats[kDuplicates];
	if ( al.IsFailedQC()  ) ++m_stats[kFailedQC];
	if ( al.IsMapped()    ) ++m_stats[kMappedReads];

	// increment strand counters
	if ( al.IsReverseStrand() ) 
		++m_stats[kReverseStrands];
	else 
		++m_stats[kForwardStrands];

	// if alignment is paired-end
	if ( al.IsPaired() ) {

		// increment PE counter
		++m_stats[kPairedEndReads];

		// increment first mate/second mate counters
		if ( al.IsFirstMate()  ) ++m_stats[kFirstMates];
		if ( al.IsSecondMate() ) ++m_stats[kSecondMates];

		// if alignment is mapped, check mate status
		if ( al.IsMapped() ) {
			// if mate mapped
			if ( al.IsMateMapped() )  {
				++m_stats[kBothMatesMapped];
			}
			// else singleton
			else 
				++m_stats[kSingletons];
		}

		// check for explicit proper pair flag
		if ( al.IsProperPair() )
			++m_stats[kProperPairs];
	}
}

void BasicStatsCollector::appendJsonImpl(json_t * jsonRootObj) {
	StatMapT::iterator iter;
	for(iter = m_stats.begin(); iter != m_stats.end(); iter++) {
		json_object_set_new(jsonRootObj, iter->first.c_str(), json_integer(iter->second));
	}
}
