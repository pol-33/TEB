# ===========================================================================
#  TEB — Sequence Mapper build system
#
#  Targets:
#    make              → build indexer + mapper (release)
#    make debug        → build with -g and -DDEBUG
#    make teb          → legacy single-binary (parsers only)
#    make exercises    → standalone exercise binaries
#    make clean        → remove all build artefacts
# ===========================================================================

CXX        = g++
CXXFLAGS   = -Wall -Wextra -std=c++17 -O3 -march=native \
             -I src/index -I src/genome -I src/io -I src/util -I vendor
DEBUGFLAGS = -DDEBUG -g

# ---------- build directory ------------------------------------------------ #
BUILDDIR = build

# ---------- shared core library objects ------------------------------------ #
# Index, genome, and I/O layers — linked into both indexer and mapper.
CORE_SRC = src/index/index.cpp \
           src/index/kmer_index.cpp \
           src/genome/genome.cpp \
           src/io/fasta.cpp \
           src/io/fastq.cpp \
           src/io/sam.cpp

# Map each src/subdir/file.cpp → build/subdir/file.o
CORE_OBJ = $(patsubst src/%.cpp,$(BUILDDIR)/%.o,$(CORE_SRC))

# ---------- per-target sources --------------------------------------------- #
INDEXER_SRC = src/indexer.cpp
MAPPER_SRC  = src/mapper.cpp
INDEXER_OBJ = $(BUILDDIR)/indexer.o
MAPPER_OBJ  = $(BUILDDIR)/mapper.o

INDEXER = $(BUILDDIR)/indexer
MAPPER  = $(BUILDDIR)/mapper

# Legacy single binary (parsers-only tool)
TEB_SRC = src/teb.cpp
TEB_OBJ = $(BUILDDIR)/teb.o
TEB     = teb.exe

# ===========================================================================
#  Top-level targets
# ===========================================================================

.PHONY: all release debug teb clean re exercises

# Default: build both new executables
all: release

release: $(INDEXER) $(MAPPER)

debug: CXXFLAGS += $(DEBUGFLAGS)
debug: $(INDEXER) $(MAPPER)

# Legacy target
teb: $(TEB)

# ---------- link rules ----------------------------------------------------- #

$(INDEXER): $(CORE_OBJ) $(INDEXER_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(MAPPER): $(CORE_OBJ) $(MAPPER_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(TEB): $(CORE_OBJ) $(TEB_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

# ---------- generic compilation rule --------------------------------------- #
# Handles both src/file.cpp → build/file.o
# and      src/subdir/file.cpp → build/subdir/file.o

$(BUILDDIR)/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create build directory if missing
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# ---------- clean ---------------------------------------------------------- #

clean:
	rm -rf $(BUILDDIR) $(TEB)

re: clean all

# ---------- standalone exercises ------------------------------------------- #

exercises: exercises/ex2_1 exercises/ex5_2

exercises/ex2_1: exercises/ex2_1.exact_search.cc
	$(CXX) $(CXXFLAGS) -o $@ $<

exercises/ex5_2: exercises/ex5_2.full_boyer_moore.cc
	$(CXX) $(CXXFLAGS) -o $@ $<