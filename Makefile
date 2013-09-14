BAMTOOLS=dep/bamtools

CFLAGS=-I$(BAMTOOLS)/src -I$(BAMTOOLS) -g
LDFLAGS=-L$(BAMTOOLS)/lib -lbamtools

bamstatsAlive: main.cc
	$(CXX) $(CFLAGS) -o bamstatsAlive main.cc $(BAMTOOLS)/src/utils/bamtools_pileup_engine.cpp $(LDFLAGS)

test1:
	curl --silent "http://bammerger.iobio.io/?binary=true&cmd=11:10108473-10188473%20http://s3.amazonaws.com/1000genomes/data/NA06984/alignment/NA06984.chrom11.ILLUMINA.bwa.CEU.low_coverage.20111114.bam%20http://s3.amazonaws.com/1000genomes/data/NA06985/alignment/NA06985.chrom11.ILLUMINA.bwa.CEU.low_coverage.20111114.bam" | ./bamstatsAlive

test2:
	curl --silent "http://bammerger.iobio.io/?binary=true&cmd=11:10108473-10188473%20http://s3.amazonaws.com/1000genomes/data/NA06984/alignment/NA06984.chrom11.ILLUMINA.bwa.CEU.low_coverage.20111114.bam%20http://s3.amazonaws.com/1000genomes/data/NA06985/alignment/NA06985.chrom11.ILLUMINA.bwa.CEU.low_coverage.20111114.bam" | ./bamstatsAlive -u 100

test3:
	curl --silent "http://bammerger.iobio.io/?binary=true&cmd=11:10108473-10188473%20http://s3.amazonaws.com/1000genomes/data/NA06984/alignment/NA06984.chrom11.ILLUMINA.bwa.CEU.low_coverage.20111114.bam%20http://s3.amazonaws.com/1000genomes/data/NA06985/alignment/NA06985.chrom11.ILLUMINA.bwa.CEU.low_coverage.20111114.bam" | ./bamstatsAlive -s 10108473 -l 80000
