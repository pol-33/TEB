CXX      = g++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra -MMD -MP -Iinclude
LDFLAGS  =

SRC_DIR  = src
OBJ_DIR  = obj

COMMON_SRCS = seed_index.cpp alignment.cpp fastq_reader.cpp
COMMON_OBJS = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(COMMON_SRCS))
DEPS = $(COMMON_OBJS:.o=.d) $(OBJ_DIR)/indexer.d $(OBJ_DIR)/mapper.d

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

-include $(DEPS)
