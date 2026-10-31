# LinearNupack DP03 binary (unity build: cpp/main.cpp #includes every part).
# Flags of the paper-v1 build. No -ffast-math / FMA contraction: results must match bit for bit.
CXX      ?= g++
CXXFLAGS ?= -O3 -std=c++17 -ffp-contract=off -Wall -Wextra
BIN      ?= build/dp03
# NUPACK rna1999 parameters (tracked).
PARAM_FILE ?= dp03/params/rna1999.dG

SRCS := $(shell find cpp -name '*.cpp' -o -name '*.hpp')

all: $(BIN)

$(BIN): $(SRCS)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -o $@ cpp/main.cpp

# V and G outside identity at the paper preset, b_V and b_G uncapped (about 1 min, 0.6 GB).
check: $(BIN)
	tests/identity_check.sh $(BIN) $(PARAM_FILE)

clean:
	rm -f $(BIN)

.PHONY: all check clean
