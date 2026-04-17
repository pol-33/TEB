CXX      = g++
UNAME_S  := $(shell uname -s)
UNAME_M  := $(shell uname -m)

BASE_CXXFLAGS = -O3 -std=c++17 -Wall -Wextra -Wpedantic -MMD -MP -Iinclude
PORTABLE_CXXFLAGS = $(BASE_CXXFLAGS)
NATIVE_CXXFLAGS   = $(BASE_CXXFLAGS)

ifeq ($(UNAME_S),Darwin)
  NATIVE_CXXFLAGS += -mcpu=apple-m1
else
  NATIVE_CXXFLAGS += -march=native
endif

CXXFLAGS ?= $(NATIVE_CXXFLAGS)

SRC_DIR = src
OBJ_DIR = obj

COMMON_SRCS = alignment.cpp buffered_io.cpp index.cpp reference.cpp search.cpp verifier.cpp
COMMON_OBJS = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(COMMON_SRCS))
DEPS        = $(COMMON_OBJS:.o=.d) $(OBJ_DIR)/indexer.d $(OBJ_DIR)/mapper.d

all: indexer mapper

portable:
	$(MAKE) CXXFLAGS='$(PORTABLE_CXXFLAGS)' all

obj:
	mkdir -p $(OBJ_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | obj
	$(CXX) $(CXXFLAGS) -c $< -o $@

indexer: $(OBJ_DIR)/indexer.o $(COMMON_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

mapper: $(OBJ_DIR)/mapper.o $(COMMON_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

test: all
	./indexer -R test/smoke_genome.fa -I test/smoke.idx
	./mapper -I test/smoke.idx -i test/smoke_reads.fastq -o test/smoke.sam -k 1
	cat test/smoke.sam

clean:
	rm -rf $(OBJ_DIR) indexer mapper test/smoke.idx test/smoke.sam

.PHONY: all portable test clean

-include $(DEPS)
