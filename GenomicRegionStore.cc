#include "GenomicRegionStore.h"

using namespace BamstatsAlive;

GenomicRegionStore::GenomicRegionStore(const std::string& regionJson) {

	json_t *root;
	json_error_t error;

	root = json_loads(regionJson.c_str(), 0, &error);

	if(!root) throw new InvalidJsonStringException;
	if(!json_is_array(root)) {
		json_decref(root);
		throw new JsonRootNotArrayException;
	}

	for(int i=0; i<json_array_size(root); i++) {
		json_t *arrayItem, *jvChr, *jvStartPos, *jvEndPos;

		arrayItem = json_array_get(root, i);
		if(!json_is_object(arrayItem)) {
			json_decref(root);
			throw new ArrayItemsNotObjectException;
		}

		jvChr = json_object_get(arrayItem, "chr");
		jvStartPos = json_object_get(arrayItem, "start");
		jvEndPos = json_object_get(arrayItem, "end");

		if(!json_is_string(jvChr) || !json_is_integer(jvStartPos) || !json_is_integer(jvEndPos)) {
			json_decref(root);
			throw new UnexpectedFieldDataTypeException;
		}

		const char *chr = json_string_value(jvChr);
		int32_t startPos = json_integer_value(jvStartPos);
		int32_t endPos = json_integer_value(jvEndPos);

		struct _regionT newRegion(chr, startPos, endPos);
		_regions.push_back(newRegion);
	}

	json_decref(root);
}

const GenomicRegionStore::GenomicRegionT& GenomicRegionStore::kRegionNotFound() {
	static GenomicRegionT notfound("", 0, 0);
	return notfound;
}

const GenomicRegionStore::GenomicRegionT& GenomicRegionStore::locateRegion(const char *chrom, int32_t pos) {

	std::vector<GenomicRegionT>::iterator it;
	for(it = _regions.begin(); it != _regions.end(); it++) {
		if(it->contains(chrom, pos))
			return *it;
	}

	return GenomicRegionStore::kRegionNotFound();
}
