CXX      = g++
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra -MMD -MP -Iinclude -march=native
LDFLAGS  =

# Enable OpenMP for parallel processing (default: enabled)
# macOS: requires llvm/clang with OpenMP (brew install libomp)
# Disable with: make OPENMP=0
OPENMP ?= 1
ifeq ($(OPENMP),1)
  # Detect macOS vs Linux
  UNAME_S := $(shell uname -s)
  ifeq ($(UNAME_S),Darwin)
    # macOS with Homebrew libomp
    CXXFLAGS += -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include
    LDFLAGS  += -L/opt/homebrew/opt/libomp/lib -lomp
  else
    # Linux with GCC
    CXXFLAGS += -fopenmp
    LDFLAGS  += -fopenmp
  endif
endif

# Enable DEBUG mode with: make DEBUG=1
ifdef DEBUG
CXXFLAGS += -DDEBUG -g
endif

SRC_DIR  = src
OBJ_DIR  = obj

COMMON_SRCS = bwt.cpp fm_index.cpp fm_search.cpp alignment.cpp fasta_reader.cpp fastq_reader.cpp
COMMON_OBJS = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(COMMON_SRCS))
DEPS = $(COMMON_OBJS:.o=.d) $(OBJ_DIR)/indexer.d $(OBJ_DIR)/mapper.d

all: indexer mapper

obj:
	mkdir -p $(OBJ_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | obj
	$(CXX) $(CXXFLAGS) -c $< -o $@

indexer: $(OBJ_DIR)/indexer.o $(COMMON_OBJS)
	$(CXX) $(LDFLAGS) $^ -o $@

mapper: $(OBJ_DIR)/mapper.o $(COMMON_OBJS)
	$(CXX) $(LDFLAGS) $^ -o $@

test: all
	./indexer -R test/tiny_genome.fa -I test/tiny.idx 2>&1
	./mapper -I test/tiny.idx -i test/tiny_reads.fastq -o test/out.sam -k 1
	cat test/out.sam

clean:
	rm -rf $(OBJ_DIR) indexer mapper test/tiny.idx test/out.sam

.PHONY: all clean test

-include $(DEPS)
