# Variables
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17 -O3
TARGET = teb
SRC = fasta_parser.cpp teb.cpp
OBJ = $(SRC:.cpp=.o)

# Regla per defecte
all: $(TARGET)

# Enllaçar l'executable
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJ)

# Compilar el fitxer .cpp a .o
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Netejar
clean:
	rm -f $(OBJ) $(TARGET)

re: clean all