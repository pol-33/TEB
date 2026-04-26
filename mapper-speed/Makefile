CXX      = g++
UNAME_S  := $(shell uname -s)
UNAME_M  := $(shell uname -m)

BASE_CXXFLAGS = -O3 -DNDEBUG -std=c++17 -Wall -Wextra -Wpedantic -MMD -MP -Iinclude \
                -fomit-frame-pointer -fstrict-aliasing
BASE_LDFLAGS  =
PORTABLE_CXXFLAGS = $(BASE_CXXFLAGS)
PORTABLE_LDFLAGS  = $(BASE_LDFLAGS)
NATIVE_CXXFLAGS   = $(BASE_CXXFLAGS)
NATIVE_LDFLAGS    = $(BASE_LDFLAGS)

ifeq ($(UNAME_S),Darwin)
  NATIVE_CXXFLAGS += -mcpu=apple-m1
else
  NATIVE_CXXFLAGS += -march=native -mtune=native -flto -fno-semantic-interposition
  NATIVE_LDFLAGS  += -flto
endif

CXXFLAGS ?= $(NATIVE_CXXFLAGS)
LDFLAGS  ?= $(NATIVE_LDFLAGS)

SRC_DIR = src
OBJ_DIR = obj

COMMON_SRCS = alignment.cpp buffered_io.cpp index.cpp reference.cpp search.cpp verifier.cpp
COMMON_OBJS = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(COMMON_SRCS))
DEPS        = $(COMMON_OBJS:.o=.d) $(OBJ_DIR)/indexer.d $(OBJ_DIR)/mapper.d

all: indexer mapper

portable:
	$(MAKE) CXXFLAGS='$(PORTABLE_CXXFLAGS)' LDFLAGS='$(PORTABLE_LDFLAGS)' all

pgo-generate:
	$(MAKE) clean
	$(MAKE) CXXFLAGS='$(NATIVE_CXXFLAGS) -fprofile-generate' LDFLAGS='$(NATIVE_LDFLAGS) -fprofile-generate' all

pgo-use:
	$(MAKE) CXXFLAGS='$(NATIVE_CXXFLAGS) -fprofile-use -fprofile-correction' \
	        LDFLAGS='$(NATIVE_LDFLAGS) -fprofile-use -fprofile-correction' all

obj:
	mkdir -p $(OBJ_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | obj
	$(CXX) $(CXXFLAGS) -c $< -o $@

indexer: $(OBJ_DIR)/indexer.o $(COMMON_OBJS)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

mapper: $(OBJ_DIR)/mapper.o $(COMMON_OBJS)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

test: all
	./indexer -R test/smoke_genome.fa -I test/smoke.idx
	./mapper -I test/smoke.idx -i test/smoke_reads.fastq -o test/smoke.sam -k 1
	cat test/smoke.sam

clean:
	rm -rf $(OBJ_DIR) indexer mapper test/smoke.idx test/smoke.sam

.PHONY: all portable pgo-generate pgo-use test clean

-include $(DEPS)
