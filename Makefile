CFLAGS=-I$(BAMTOOLS)/src -I$(BAMTOOLS) -Ilib/jansson-2.5/src -g
LDFLAGS=-L$(BAMTOOLS)/lib -lbamtools

SOURCES=main.cc \
		$(BAMTOOLS)/src/utils/bamtools_pileup_engine.cpp

STATLIBS=lib/jansson-2.5/src/.libs/libjansson.a

bamstatsAlive: checkvar libjansson main.cc
	$(CXX) $(CFLAGS) -o bamstatsAlive $(SOURCES) $(STATLIBS) $(LDFLAGS)

checkvar:
	@if [ "x$(BAMTOOLS)" = "x" ]; then echo "BAMTOOLS need to be defined"; exit 1; fi

libjansson: lib/jansson-2.5/src/.libs/libjansson.a

lib/jansson-2.5/src/.libs/libjansson.a:
	@if [ ! -d lib ]; then mkdir lib; fi
	@if [ ! -d lib/jansson-2.5 ]; then cd lib; curl -o - http://www.digip.org/jansson/releases/jansson-2.5.tar.gz | tar -xzf - ; fi
	@cd lib/jansson-2.5; ./configure --disable-shared --enable-static; make; cd ../..

test1:
	curl --silent "http://bammerger.iobio.io/?binary=true&cmd=11:10108473-10188473%20http://s3.amazonaws.com/1000genomes/data/NA06984/alignment/NA06984.chrom11.ILLUMINA.bwa.CEU.low_coverage.20111114.bam%20http://s3.amazonaws.com/1000genomes/data/NA06985/alignment/NA06985.chrom11.ILLUMINA.bwa.CEU.low_coverage.20111114.bam" | ./bamstatsAlive

test2:
	curl --silent "http://bammerger.iobio.io/?binary=true&cmd=11:10108473-10188473%20http://s3.amazonaws.com/1000genomes/data/NA06984/alignment/NA06984.chrom11.ILLUMINA.bwa.CEU.low_coverage.20111114.bam%20http://s3.amazonaws.com/1000genomes/data/NA06985/alignment/NA06985.chrom11.ILLUMINA.bwa.CEU.low_coverage.20111114.bam" | ./bamstatsAlive -u 100

test3:
	curl --silent "http://bammerger.iobio.io/?binary=true&cmd=11:10108473-10188473%20http://s3.amazonaws.com/1000genomes/data/NA06984/alignment/NA06984.chrom11.ILLUMINA.bwa.CEU.low_coverage.20111114.bam%20http://s3.amazonaws.com/1000genomes/data/NA06985/alignment/NA06985.chrom11.ILLUMINA.bwa.CEU.low_coverage.20111114.bam" | ./bamstatsAlive -s 10108473 -l 80000
