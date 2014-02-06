CXXFLAGS=-I$(BAMTOOLS)/src -I$(BAMTOOLS) -I../lib/jansson-2.5/src -I..
LDFLAGS=-L$(BAMTOOLS)/lib -lbamtools

.SUFFIXES: .cc

TEST_SOURCES=testGenomicRegionStore.cc

TEST_OBJECTS=$(TEST_SOURCES:.cc=.o)

TESTS=$(TEST_OBJECTS:.o=)

STATLIBS=../lib/jansson-2.5/src/.libs/libjansson.a

PARENT_OBJECTS=$(wildcard ../*.o)
LIB_PARENT_OBJECTS=$(filter-out ../main.o, $(PARENT_OBJECTS))

.PHONY: all clean test

all: test

clean:
	rm -rf $(TEST_OBJECTS) $(TEST_BINS)

.cc.o:
	$(CXX) -c $< $(CXXFLAGS)

$(TESTS) : $(TEST_OBJECTS) $(LIB_PARENT_OBJECTS)
	$(CXX) -o $@ $< $(LIB_PARENT_OBJECTS) $(STATLIBS) $(LDFLAGS)

test: $(TESTS)
	@./run-test.sh $<
