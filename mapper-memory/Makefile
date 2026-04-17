CXX      = g++
BASE_CXXFLAGS = -O3 -std=c++17 -Wall -Wextra -MMD -MP -Iinclude
CXXFLAGS = $(BASE_CXXFLAGS)
LDFLAGS  =

SIMD ?= auto
ARCH := $(shell uname -m)
AVX512_FLAG_CHECK := $(shell printf 'int main(void){return 0;}\n' | $(CXX) -x c++ -std=c++17 -mavx512f -mavx512bw -c -o /dev/null - >/dev/null 2>&1 && echo 1 || echo 0)
HAVE_AVX512_IMPL := 0

ifeq ($(SIMD),auto)
  ifeq ($(ARCH),x86_64)
    ifeq ($(AVX512_FLAG_CHECK),1)
      HAVE_AVX512_IMPL := 1
    endif
  endif
else ifeq ($(SIMD),off)
  HAVE_AVX512_IMPL := 0
else ifeq ($(SIMD),avx512)
  ifneq ($(ARCH),x86_64)
    $(error SIMD=avx512 requires an x86_64 build host)
  endif
  ifneq ($(AVX512_FLAG_CHECK),1)
    $(error SIMD=avx512 requires compiler support for -mavx512f -mavx512bw)
  endif
  HAVE_AVX512_IMPL := 1
else
  $(error SIMD must be one of: auto, off, avx512)
endif

CXXFLAGS += -DMAPPER_MEMORY_HAVE_AVX512_IMPL=$(HAVE_AVX512_IMPL)
AVX512_CXXFLAGS = -mavx512f -mavx512bw

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

COMMON_SRCS = bwt.cpp fm_index.cpp fm_search.cpp alignment.cpp fasta_reader.cpp fastq_reader.cpp simd_dispatch.cpp
COMMON_OBJS = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(COMMON_SRCS))
SIMD_AVX512_OBJS =
ifeq ($(HAVE_AVX512_IMPL),1)
  SIMD_AVX512_OBJS += $(OBJ_DIR)/simd_avx512.o
endif
DEPS = $(COMMON_OBJS:.o=.d) $(SIMD_AVX512_OBJS:.o=.d) $(OBJ_DIR)/indexer.d $(OBJ_DIR)/mapper.d $(OBJ_DIR)/simd_selftest.d

all: indexer mapper

obj:
	mkdir -p $(OBJ_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | obj
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/simd_avx512.o: $(SRC_DIR)/simd_avx512.cpp | obj
	$(CXX) $(CXXFLAGS) $(AVX512_CXXFLAGS) -c $< -o $@

indexer: $(OBJ_DIR)/indexer.o $(COMMON_OBJS) $(SIMD_AVX512_OBJS)
	$(CXX) $(LDFLAGS) $^ -o $@

mapper: $(OBJ_DIR)/mapper.o $(COMMON_OBJS) $(SIMD_AVX512_OBJS)
	$(CXX) $(LDFLAGS) $^ -o $@

simd_selftest: $(OBJ_DIR)/simd_selftest.o $(COMMON_OBJS) $(SIMD_AVX512_OBJS)
	$(CXX) $(LDFLAGS) $^ -o $@

test: all simd_selftest
	MAPPER_SIMD=off ./simd_selftest
	MAPPER_SIMD=auto ./simd_selftest
	./indexer -R test/tiny_genome.fa -I test/tiny.idx
	MAPPER_SIMD=off ./mapper -I test/tiny.idx -i test/tiny_reads.fastq -o test/out.scalar.sam -k 1
	cmp -s test/out.scalar.sam test/out_L1.sam
	MAPPER_SIMD=auto ./mapper -I test/tiny.idx -i test/tiny_reads.fastq -o test/out.auto.sam -k 1
	cmp -s test/out.auto.sam test/out_L1.sam

clean:
	rm -rf $(OBJ_DIR) indexer mapper simd_selftest test/tiny.idx test/out.sam test/out.scalar.sam test/out.auto.sam

.PHONY: all clean test simd_selftest

-include $(DEPS)
