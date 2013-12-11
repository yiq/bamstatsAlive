#ifndef ABSTRACTSTATCOLLECTOR_H
#define ABSTRACTSTATCOLLECTOR_H

#pragma once

#include <vector>
#include <api/BamAlignment.h>
#include <jansson.h>

namespace BamstatsAlive {

	class AbstractStatCollector;

	typedef std::vector<const AbstractStatCollector *> StatCollectorPtrVec;

	class AbstractStatCollector {
		protected:
			StatCollectorPtrVec _children;

		public:
			AbstractStatCollector();
			virtual ~AbstractStatCollector();

			void addChild(const AbstractStatCollector * child);
			void removeChild(const AbstractStatCollector * child);

			virtual void processAlignment(const BamTools::BamAlignment& al, const BamTools::RefVector& refVector) = 0;
			virtual json_t * appendJson(json_t * jsonRootObj) = 0;
	};

}

#endif
