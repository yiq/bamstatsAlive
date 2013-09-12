BAMTOOLS=dep/bamtools

CFLAGS=-I$(BAMTOOLS)/include
LDFLAGS=-L$(BAMTOOLS)/lib -lbamtools

bamStatsAlive: main.cc
	$(CXX) $(CFLAGS) -o bamStatsAlive main.cc $(LDFLAGS)
