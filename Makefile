CFLAGS=-std=c++11 -I$(BAMTOOLS)/src -I$(BAMTOOLS) -Ilib/jansson-2.8/src 
LDFLAGS=-L$(BAMTOOLS)/lib -L$(BAMTOOLS)/build/src/api -lbamtools -lz

.SUFFIXES: .cc

SOURCES=main.cc \
		AbstractStatCollector.cc \
		BasicStatsCollector.cc \
		HistogramStatsCollector.cc \
		CoverageMapStatsCollector.cc \
		GenomicRegionStore.cc

OBJECTS=$(SOURCES:.cc=.o)
OBJECTS+=bamtools_pileup_engine.o 


#OBJECTS=main.o \
#		AbstractStatCollector.o \
#		BasicStatsCollector.o \
#		HistogramStatsCollector.o \
#		CoverageMapStatsCollector.o \
#		GenomicRegionStore.o \
#		bamtools_pileup_engine.o 

STATLIBS=lib/jansson-2.8/src/.libs/libjansson.a

all: release

debug: CFLAGS += -DDEBUG -g -pg
debug: bamstatsAlive

release: CFLAGS += -DRELEASE -O2
release: bamstatsAlive

clean:
	rm -rf *.o *.dSYM bamstatsAlive bamstatsAliveCommon.hpp.gch

clean-dep:
	make -C lib/jansson-2.8 clean

.PHONY: all clean clean-dep

.cc.o :
	$(CXX) -c $< $(CFLAGS) -include bamstatsAliveCommon.hpp

bamtools_pileup_engine.o: $(BAMTOOLS)/src/utils/bamtools_pileup_engine.cpp
	$(CXX) -c $< $(CFLAGS)

bamstatsAliveCommon.hpp.gch: bamstatsAliveCommon.hpp
	$(CXX) $(CFLAGS) -x c++-header $< -Winvalid-pch -o $@

bamstatsAlive: checkvar libjansson bamstatsAliveCommon.hpp.gch $(OBJECTS)
	$(CXX) $(CFLAGS) -o bamstatsAlive $(OBJECTS) $(STATLIBS) $(LDFLAGS)

checkvar:
	@if [ "x$(BAMTOOLS)" = "x" ]; then echo "BAMTOOLS need to be defined"; exit 1; fi

libjansson: lib/jansson-2.8/src/.libs/libjansson.a

lib/jansson-2.8/src/.libs/libjansson.a:
	@if [ ! -d lib ]; then mkdir lib; fi
	@if [ ! -d lib/jansson-2.8 ]; then cd lib; curl -Lo - https://digip.org/jansson/releases/jansson-2.8.tar.gz | tar -xzf - ; fi
	@cd lib/jansson-2.8; ./configure --disable-shared --enable-static; make; cd ../..

test1:
	curl --silent "http://bammerger.iobio.io/?binary=true&cmd=11:10108473-10188473%20http://s3.amazonaws.com/1000genomes/data/NA06984/alignment/NA06984.chrom11.ILLUMINA.bwa.CEU.low_coverage.20111114.bam%20http://s3.amazonaws.com/1000genomes/data/NA06985/alignment/NA06985.chrom11.ILLUMINA.bwa.CEU.low_coverage.20111114.bam" | ./bamstatsAlive

test2:
	curl --silent "http://bammerger.iobio.io/?binary=true&cmd=11:10108473-10188473%20http://s3.amazonaws.com/1000genomes/data/NA06984/alignment/NA06984.chrom11.ILLUMINA.bwa.CEU.low_coverage.20111114.bam%20http://s3.amazonaws.com/1000genomes/data/NA06985/alignment/NA06985.chrom11.ILLUMINA.bwa.CEU.low_coverage.20111114.bam" | ./bamstatsAlive -u 100

test3:
	curl --silent "http://bammerger.iobio.io/?binary=true&cmd=11:10108473-10188473%20http://s3.amazonaws.com/1000genomes/data/NA06984/alignment/NA06984.chrom11.ILLUMINA.bwa.CEU.low_coverage.20111114.bam%20http://s3.amazonaws.com/1000genomes/data/NA06985/alignment/NA06985.chrom11.ILLUMINA.bwa.CEU.low_coverage.20111114.bam" | ./bamstatsAlive -s 10108473 -l 80000
