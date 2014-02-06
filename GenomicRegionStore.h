#ifndef GENOMICREGIONSTORE_H
#define GENOMICREGIONSTORE_H

#pragma once

#include <stdint.h>
#include <vector>
#include <iostream>

namespace BamstatsAlive {


	class GenomicRegionStore {

		public:
			typedef struct _regionT {
				const char * chrom;
				int32_t startPos;
				int32_t endPos;

				_regionT(const char *chrom, int32_t startPos, int32_t endPos) : 
					chrom(NULL), startPos(startPos), endPos(endPos) 
				{ this->chrom = strdup(chrom); }

				bool contains(const char *chrom, int32_t pos) const {
					if(strcmp(chrom, this->chrom) != 0) return false;
					if(pos < startPos || pos > endPos) return false;
					return true;
				}
			} GenomicRegionT;

		protected:
			std::vector<GenomicRegionT> _regions;

		public:
			GenomicRegionStore(const std::string& regionJson);

			inline const std::vector<GenomicRegionT>& regions() { return _regions; }

			// methods for locating a region
			static const GenomicRegionT& kRegionNotFound(); 
			const GenomicRegionT& locateRegion(const char *chrom, int32_t pos);

			
			class InvalidJsonStringException {};
			class JsonRootNotArrayException {};
			class ArrayItemsNotObjectException {};
			class UnexpectedFieldDataTypeException {};
	};
}

#endif
