# Variables
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17 -O3 -I include/
TARGET = teb
SRC = src/fasta_parser.cpp src/fastq_parser.cpp src/teb.cpp
OBJ = $(SRC:src/%.cpp=src/%.o)

# Regla per defecte
all: $(TARGET)

# Enllaçar l'executable
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJ)

# Compilar el fitxer .cpp a .o
src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Netejar
clean:
	rm -f $(OBJ) $(TARGET)

re: clean all

# Standalone exercises (compiled separately)
exercises/ex2_1: exercises/ex2_1.exact_search.cc
	$(CXX) $(CXXFLAGS) -o $@ $<

exercises/ex5_2: exercises/ex5_2.full_boyer_moore.cc
	$(CXX) $(CXXFLAGS) -o $@ $<