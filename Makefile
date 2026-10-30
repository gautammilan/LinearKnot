# LinearNupack DP03 binary (unity build: cpp/main.cpp #includes every part).
# Flags of the paper-v1 build. No -ffast-math / FMA contraction: results must match bit for bit.
CXX      ?= g++
CXXFLAGS ?= -O3 -std=c++17 -ffp-contract=off -Wall -Wextra
BIN      ?= build/dp03

SRCS := $(shell find cpp -name '*.cpp' -o -name '*.hpp')

all: $(BIN)

$(BIN): $(SRCS)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -o $@ cpp/main.cpp

clean:
	rm -f $(BIN)

.PHONY: all clean
