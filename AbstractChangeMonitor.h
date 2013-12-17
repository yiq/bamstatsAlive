#ifndef ABSTRACTCHANGEMONITOR_H
#define ABSTRACTCHANGEMONITOR_H

#pragma once

namespace BamstatsAlive {

	template<class T>
	class AbstractChangeMonitor {
		public:
			virtual void addValue(T value) = 0;
			virtual bool isSatisfied() = 0;
	};
}

#endif
