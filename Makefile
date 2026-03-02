# Variables
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17 -O3
DEBUGFLAGS = -DDEBUG -g

TARGET = teb.exe
SRC = fasta_parser.cpp teb.cpp utils.cpp
OBJ = $(SRC:.cpp=.o)

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
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean
clean:
	rm -f $(OBJ) $(TARGET)

re: clean all