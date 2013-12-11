#include "AbstractStatCollector.h"

#include <vector>
#include <algorithm>

using namespace BamstatsAlive;

AbstractStatCollector::AbstractStatCollector() {

}

AbstractStatCollector::~AbstractStatCollector() {

}

void AbstractStatCollector::addChild(const AbstractStatCollector *child) {

	// Make sure that the input is good
	if(child == NULL) return;

	// Make sure the input is not yet a child
	StatCollectorPtrVec::const_iterator loc = std::find(_children.begin(), _children.end(), child);
	if(loc != _children.end()) return;

	// Insert the input to the end of the children list
	_children.push_back(child);
}

void AbstractStatCollector::removeChild(const AbstractStatCollector * child) {

	// Make sure that the input is good
	if(child == NULL) return;

	// Make sure the input is a child
	StatCollectorPtrVec::const_iterator loc = std::find(_children.begin(), _children.end(), child);
	if(loc == _children.end()) return;

	_children.erase(loc);
}
