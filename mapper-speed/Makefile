CXX      = g++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra -Iinclude

all: indexer mapper

indexer: src/indexer.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

mapper: src/mapper.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f indexer mapper

.PHONY: all clean
