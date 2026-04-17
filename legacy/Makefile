# Variables
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17 -O3 -I include/
DEBUGFLAGS = -DDEBUG -g

TARGET = teb.exe
SRC = src/fasta_parser.cpp src/fastq_parser.cpp src/utils.cpp src/teb.cpp
OBJ = $(SRC:src/%.cpp=src/%.o)

# Default target = release
all: release

# Release build
release: $(TARGET)

# Debug build
debug: CXXFLAGS += $(DEBUGFLAGS)
debug: $(TARGET)

# Link executable
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJ)

# Compile .cpp to .o
src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean
clean:
	rm -f $(OBJ) $(TARGET)

re: clean all

# Standalone exercises (compiled separately)
exercises/ex2_1: exercises/ex2_1.exact_search.cc
	$(CXX) $(CXXFLAGS) -o $@ $<

exercises/ex5_2: exercises/ex5_2.full_boyer_moore.cc
	$(CXX) $(CXXFLAGS) -o $@ $<
