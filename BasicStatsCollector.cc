#include "BasicStatsCollector.h"
#include "StandardDeviationChangeMonitor.h"
#include "DeltaAverageRatioChangeMonitor.h"
#include <iostream>

using namespace BamstatsAlive;

BasicStatsCollector::BasicStatsCollector() {
	_stats[kTotalReads] = 0;
	_stats[kMappedReads] = 0;
	_stats[kForwardStrands] = 0;
	_stats[kReverseStrands] = 0;
	_stats[kFailedQC] = 0;
	_stats[kDuplicates] = 0;
	_stats[kPairedEndReads] = 0;
	_stats[kProperPairs] = 0;
	_stats[kBothMatesMapped] = 0;
	_stats[kFirstMates] = 0;
	_stats[kSecondMates] = 0;
	_stats[kSingletons] = 0;

	_stats.clear();

	StatMapT::iterator iter;
	for(iter = _stats.begin(); iter != _stats.end(); iter++) {
		std::cerr<<"Initializing: "<<iter->first<<std::endl;
	}

	_changeMonitor = new StandardDeviationChangeMonitor<double>(3, 0.05);

}

void BasicStatsCollector::processAlignmentImpl(const BamTools::BamAlignment& al, const BamTools::RefVector& refVector) {
	// increment total alignment counter
	++_stats[kTotalReads];

	// incrememt counters for pairing-independent flags
	if ( al.IsDuplicate() ) ++_stats[kDuplicates];
	if ( al.IsFailedQC()  ) ++_stats[kFailedQC];
	if ( al.IsMapped()    ) ++_stats[kMappedReads];

	// increment strand counters
	if ( al.IsReverseStrand() ) 
		++_stats[kReverseStrands];
	else 
		++_stats[kForwardStrands];

	// if alignment is paired-end
	if ( al.IsPaired() ) {

		// increment PE counter
		++_stats[kPairedEndReads];

		// increment first mate/second mate counters
		if ( al.IsFirstMate()  ) ++_stats[kFirstMates];
		if ( al.IsSecondMate() ) ++_stats[kSecondMates];

		// if alignment is mapped, check mate status
		if ( al.IsMapped() ) {
			// if mate mapped
			if ( al.IsMateMapped() )  {
				++_stats[kBothMatesMapped];
			}
			// else singleton
			else 
				++_stats[kSingletons];
		}

		// check for explicit proper pair flag
		if ( al.IsProperPair() )
			++_stats[kProperPairs];
	}
}

void BasicStatsCollector::appendJsonImpl(json_t * jsonRootObj) {
	StatMapT::iterator iter;
	for(iter = _stats.begin(); iter != _stats.end(); iter++) {
		json_object_set_new(jsonRootObj, iter->first.c_str(), json_integer(iter->second));
	}

	_changeMonitor->addValue(_stats[kMappedReads] / static_cast<double>(_stats[kTotalReads]));
	std::cerr<<"Fraction: "<<_stats[kMappedReads] / static_cast<double>(_stats[kTotalReads])<<std::endl;
	std::cerr<<"Stdev: "<<dynamic_cast<StandardDeviationChangeMonitor<double>*>(_changeMonitor)->getStdev()<<std::endl;
}
