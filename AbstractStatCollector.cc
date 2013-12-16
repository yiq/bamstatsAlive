#include "AbstractStatCollector.h"

#include <vector>
#include <algorithm>

using namespace BamstatsAlive;

AbstractStatCollector::AbstractStatCollector() {
	_children.clear();
}

AbstractStatCollector::~AbstractStatCollector() {

}

void AbstractStatCollector::addChild(AbstractStatCollector *child) {

	// Make sure that the input is good
	if(child == NULL) return;

	// Make sure the input is not yet a child
	StatCollectorPtrVec::const_iterator loc = std::find(_children.begin(), _children.end(), child);
	if(loc != _children.end()) return;

	// Insert the input to the end of the children list
	_children.push_back(child);
}

void AbstractStatCollector::removeChild(AbstractStatCollector * child) {

	// Make sure that the input is good
	if(child == NULL) return;

	// Make sure the input is a child
	StatCollectorPtrVec::const_iterator loc = std::find(_children.begin(), _children.end(), child);
	if(loc == _children.end()) return;

	_children.erase(loc);
}

void AbstractStatCollector::processAlignment(const BamTools::BamAlignment& al, const BamTools::RefVector& refVector) {
	this->processAlignmentImpl(al, refVector);

	StatCollectorPtrVec::iterator iter;
	for(iter = _children.begin(); iter != _children.end(); iter++) {
		(*iter)->processAlignment(al, refVector);
	}
}

json_t * AbstractStatCollector::appendJson(json_t * jsonRootObj) {
	if(jsonRootObj == NULL)
		jsonRootObj = json_object();

	this->appendJsonImpl(jsonRootObj);
	
	StatCollectorPtrVec::iterator iter;
	for(iter = _children.begin(); iter != _children.end(); iter++) {
		(*iter)->appendJson(jsonRootObj);
	}

	return jsonRootObj;
}
