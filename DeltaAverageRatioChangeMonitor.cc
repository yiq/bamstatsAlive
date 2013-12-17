#include "DeltaAverageRatioChangeMonitor.h"

using namespace BamstatsAlive;

template<class T>
DeltaAverageRatioChangeMonitor<T>::DeltaAverageRatioChangeMonitor(double threshold) :
	_threshold(threshold), count(0), 
	_total(0), _last(0), _delta(0) {
	
}

template<class T>
void DeltaAverageRatioChangeMonitor<T>::addValue(T value) {
	_total += value;
	_delta = value - _last;
	_last = value;
	count++;
}

template<class T>
bool DeltaAverageRatioChangeMonitor<T>::isSatisfied() {
	double average = _total / static_cast<double>(count);
	double daRatio = _delta / average;
	return daRatio < _threshold;
}
