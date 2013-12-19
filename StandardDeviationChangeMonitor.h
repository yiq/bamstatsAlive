#ifndef STANDARDDEVIATIONCHANGEMONITOR_H
#define STANDARDDEVIATIONCHANGEMONITOR_H

#pragma once

#include "AbstractChangeMonitor.h"
#include <cmath>

namespace BamstatsAlive {

	template<class T>
		class StandardDeviationChangeMonitor : public AbstractChangeMonitor<T> {
			protected:
				T * _trailingVal;
				double _threshold;
				unsigned int _count;
				unsigned int _trailLength;

				double _stdev() {
					T total = 0;
					unsigned int n = _count < _trailLength ? _count : _trailLength;
					for(int i=0; i<n; i++) 
						total += _trailingVal[i];

					double average = static_cast<double>(total) / static_cast<double>(n);
					double sqrSum = 0;
					for(int i=0; i<n; i++) 
						sqrSum += pow(static_cast<double>(_trailingVal[i]) - average, 2);

					return sqrt(sqrSum / n);
				}

			public:
				StandardDeviationChangeMonitor(unsigned int trailLength, double threshold) :
					_count(0), _trailLength(trailLength), _threshold(threshold) {
						_trailingVal = new T[trailLength];
					}

				virtual void addValue(T value) {
					if(_count >= _trailLength) {
						for(unsigned int i=0; i<_trailLength-1; i++)
							_trailingVal[i] = _trailingVal[i+1];
						_trailingVal[_trailLength - 1] = value;
						_count++;
					}
					else {
						_trailingVal[_count++] = value;
					}
				}

				virtual bool isSatisfied() {

					if(_count < _trailLength) return false;

					double stdev = _stdev();
					
					return stdev < _threshold;
				}

				double getStdev() {
					if(_count < _trailLength) return -1.0;

					return _stdev();
				}
		};
}

#endif
