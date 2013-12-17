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
			DeltaAverageRatioChangeMonitor(double threshold);
			virtual void addValue(T value);
			virtual bool isSatisfied();
	};
}

#endif
