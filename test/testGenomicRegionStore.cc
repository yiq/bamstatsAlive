#include "../GenomicRegionStore.h"

#include <string>
#include <iostream>

#define ASSERT_EQ(expr, expect, msg) { if ((expr) != (expect)) {std::cerr<<(msg)<<std::endl; exit(1);} }
#define ASSERT_UNEQ(expr, expect, msg) { if ((expr) == (expect)) {std::cerr<<(msg)<<std::endl; exit(1);} }

#define ASSERT_NOTNULL(expr, msg) ASSERT_UNEQ((expr), NULL, (msg))


using namespace std;
using namespace BamstatsAlive;

std::string testJsonStr = "[{\"start\":1,\"end\":10001,\"chr\":\"11\"},{\"start\":13500652,\"end\":13510652,\"chr\":\"11\"},{\"start\":27001304,\"end\":27011304,\"chr\":\"11\"},{\"start\":40501955,\"end\":40511955,\"chr\":\"11\"},{\"start\":54002607,\"end\":54012607,\"chr\":\"11\"},{\"start\":67503259,\"end\":67513259,\"chr\":\"11\"},{\"start\":81003910,\"end\":81013910,\"chr\":\"11\"},{\"start\":94504562,\"end\":94514562,\"chr\":\"11\"},{\"start\":108005213,\"end\":108015213,\"chr\":\"11\"},{\"start\":121505865,\"end\":121515865,\"chr\":\"11\"}]";

int main(int argc, char* argv[]) {

	GenomicRegionStore *store;

	try {
		store = new GenomicRegionStore(testJsonStr);
	}
	catch (...) {
		std::cerr<<"Exception thrown during store creation"<<std::endl;
		return 1;
	}

	ASSERT_NOTNULL(store, "Genomic Store should not be NULL");

	ASSERT_EQ(store->regions().size(), 10, "10 Genomic regions should have been observed");

	const GenomicRegionStore::GenomicRegionT& region1 = store->locateRegion("11", 500);
	ASSERT_UNEQ(&region1, &GenomicRegionStore::kRegionNotFound(), "Chrom 11 Pos 500 should be found");

	const GenomicRegionStore::GenomicRegionT& region2 = store->locateRegion("11", 20000);
	ASSERT_EQ(&region2, &GenomicRegionStore::kRegionNotFound(), "Chrom 11 Pos 20000 should not be found");

	ASSERT_EQ(region1.contains("11", 200), true, "Chrom 11 Pos 200 should be contained in the region 11:1-10001");
	ASSERT_EQ(region1.contains("11", 800), true, "Chrom 11 Pos 800 should be contained in the region 11:1-10001");
	ASSERT_EQ(region1.contains("11", 12000), false, "Chrom 11 Pos 12000 should be contained in the region 11:1-10001");

	
	const GenomicRegionStore::GenomicRegionT& readRegion1 = store->locateRegion("11", 13500579);
	ASSERT_EQ(&readRegion1, &GenomicRegionStore::kRegionNotFound(), "Start of the sample read should not be found");
	const GenomicRegionStore::GenomicRegionT& readRegion2 = store->locateRegion("11", 13500655);
	ASSERT_UNEQ(&readRegion2, &GenomicRegionStore::kRegionNotFound(), "End of the sample read should be found");


	delete store;
	return 0;
}
