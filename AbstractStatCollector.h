#ifndef ABSTRACTSTATCOLLECTOR_H
#define ABSTRACTSTATCOLLECTOR_H

#pragma once

namespace BamstatsAlive {

	class AbstractStatCollector;

	typedef std::vector<AbstractStatCollector *> StatCollectorPtrVec;

	/**
	 * The base class for all statistics collectors
	 *
	 * A statistics collector will implement two virtual functions: 
	 *   - processAlignment() to update statistics
	 *   - appendJson() to create the json representation of the statistics
	 *
	 * These statistics collectors can be organized into a tree with the
	 * addChild() and removeChild() functions. User code will only need to call
	 * the public processAlignment() and appendJson() functions on the root
	 * object, and the action will be propagated across all child nodes. The
	 * actual implementation of specific collectors is encapsulated by the
	 * protected processAlignmentImpl() and appendJsonImpl() functions
	 */
	class AbstractStatCollector {
		protected:
			StatCollectorPtrVec _children;

			/**
			 * Process the alignment and update statistics
			 *
			 * @param al The alignment read
			 * @param refVector The reference the read is aligned to
			 */
			virtual void processAlignmentImpl(const BamTools::BamAlignment& al, const BamTools::RefVector& refVector) = 0;

			/**
			 * Append statistics as json
			 *
			 * @param jsonRootObj The json root object to which the outputs are appended
			 */
			virtual void appendJsonImpl(json_t * jsonRootObj) = 0;

			/** 
			 * Check if the statistics collector is satisfied with the data it
			 * has seen so far. Note that the defualt implementation of this
			 * function always returns false, which means that by default, a
			 * collector will keep processing reads indefinitely.
			 *
			 * @return true if the data is considered to be sufficient, false otherwise
			 */
			virtual bool isSatisfiedImpl();

		public:
			AbstractStatCollector();
			virtual ~AbstractStatCollector();

			/**
			 * Add a statistics collector as the child of the current collector
			 * 
			 * @param child The child collector to be added
			 */
			void addChild(AbstractStatCollector * child);

			/**
			 * Remove a statistics collector from the children list of the current collector
			 *
			 * @param child The child collector to be removed
			 */
			void removeChild(AbstractStatCollector * child);

			/**
			 * Process an alignment by the collector tree
			 *
			 * The alignment will be passed to the processAlignmentImpl function of the
			 * current collector, and the processAlignment function of all the children
			 * collectors.
			 *
			 * @param al The alignment read
			 * @param refVector The reference the read is aligned to
			 */
			void processAlignment(const BamTools::BamAlignment& al, const BamTools::RefVector& refVector);

			/**
			 * Create json of the collector tree
			 *
			 * Json objects representing all the statistics collected by the
			 * current tree will be appended into the given json root object.
			 * If no root object is given, a new one will be created. Current
			 * collector's appendJsonImpl function and all the children's 
			 * appendJson function will be called on the root object.
			 *
			 * @param jsonRootObj The json root object to which all statistics are appended to
			 * @return The json object representing the root of all the statistics
			 */
			json_t * appendJson(json_t * jsonRootObj = NULL);

			/**
			 * Check satisfy-ness of the collector tree
			 *
			 * @return true if all collectors in the tree are satisfied, false otherwise
			 */
			bool isSatisfied();
	};

}

#endif
