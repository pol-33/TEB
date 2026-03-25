CXX      = g++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra -Iinclude
LDFLAGS  =

SRC_DIR  = src
OBJ_DIR  = obj

COMMON_SRCS = bwt.cpp fm_index.cpp backward_search.cpp alignment.cpp \
              fastq_reader.cpp fasta_reader.cpp
COMMON_OBJS = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(COMMON_SRCS))

all: indexer mapper

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

indexer: $(OBJ_DIR)/indexer.o $(COMMON_OBJS)
	$(CXX) $(LDFLAGS) $^ -o $@

mapper: $(OBJ_DIR)/mapper.o $(COMMON_OBJS)
	$(CXX) $(LDFLAGS) $^ -o $@

clean:
	rm -rf $(OBJ_DIR) indexer mapper

re: clean all

.PHONY: all clean re
