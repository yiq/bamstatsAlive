#ifndef DELTAAVERAGERATIOCHANGEMONITOR_H
#define DELTAAVERAGERATIOCHANGEMONITOR_H

#pragma once

#include "AbstractChangeMonitor.h"

namespace BamstatsAlive {

	template<class T>
	class DeltaAverageRatioChangeMonitor : public AbstractChangeMonitor<T> {
		protected:
			T _total;
			T _last;
			T _delta;
			unsigned long count;
			double _threshold;

		public:
			DeltaAverageRatioChangeMonitor(double threshold) : 
				_threshold(threshold), count(0),
				_total(0), _last(0), _delta(0) {

			}

			virtual void addValue(T value) {
				_total += value;
				_delta = value - _last;
				_last = value;
				count++;
			}

			virtual bool isSatisfied() {
				double average = _total / static_cast<double>(count);
				double daRatio = _delta / average;
				return daRatio < _threshold;
			}
	};
}

#endif
