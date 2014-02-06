#include "BasicStatsCollector.h"
#include "StandardDeviationChangeMonitor.h"
#include "DeltaAverageRatioChangeMonitor.h"

using namespace BamstatsAlive;

static const unsigned int kCMTrailLength = 5;
static const double kCMThreshold = 0.001;

BasicStatsCollector::BasicStatsCollector() {

	_stats.clear();

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
	_stats[kLastReadPos] = 0;


#ifdef DEBUG
	StatMapT::iterator iter;
	for(iter = _stats.begin(); iter != _stats.end(); iter++) {
		std::cerr<<"Initializing: "<<iter->first<<std::endl;
	}
#endif

	_monitors.clear();
	_monitors[kMappedReads] 		= new StandardDeviationChangeMonitor<double>(kCMTrailLength, kCMThreshold);
	_monitors[kForwardStrands] 		= new StandardDeviationChangeMonitor<double>(kCMTrailLength, kCMThreshold);
	_monitors[kReverseStrands] 		= new StandardDeviationChangeMonitor<double>(kCMTrailLength, kCMThreshold);
	_monitors[kFailedQC] 			= new StandardDeviationChangeMonitor<double>(kCMTrailLength, kCMThreshold);
	_monitors[kDuplicates] 			= new StandardDeviationChangeMonitor<double>(kCMTrailLength, kCMThreshold);
	_monitors[kPairedEndReads] 		= new StandardDeviationChangeMonitor<double>(kCMTrailLength, kCMThreshold);
	_monitors[kProperPairs] 		= new StandardDeviationChangeMonitor<double>(kCMTrailLength, kCMThreshold);
	_monitors[kBothMatesMapped] 	= new StandardDeviationChangeMonitor<double>(kCMTrailLength, kCMThreshold);
	_monitors[kFirstMates] 			= new StandardDeviationChangeMonitor<double>(kCMTrailLength, kCMThreshold);
	_monitors[kSecondMates] 		= new StandardDeviationChangeMonitor<double>(kCMTrailLength, kCMThreshold);
	_monitors[kSingletons] 			= new StandardDeviationChangeMonitor<double>(kCMTrailLength, kCMThreshold);
}

void BasicStatsCollector::processAlignmentImpl(const BamTools::BamAlignment& al, const BamTools::RefVector& refVector) {
	// increment total alignment counter
	++_stats[kTotalReads];

	_stats[kLastReadPos] = al.Position;

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
	StatMapT::iterator sIter;
	for(sIter = _stats.begin(); sIter != _stats.end(); sIter++) {
		json_object_set_new(jsonRootObj, sIter->first.c_str(), json_integer(sIter->second));
	}

	ChangeMonitorMapT::iterator cIter;
	for(cIter = _monitors.begin(); cIter != _monitors.end(); cIter++) {
		dynamic_cast<StandardDeviationChangeMonitor<double> *>(cIter->second)->addValue(_stats[cIter->first] / static_cast<double>(_stats[kTotalReads]));
	}

	bool consensus = true;
}
