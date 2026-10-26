CXX ?= g++
CXXFLAGS ?= -std=c++20 -O3 -Wall -Wextra -Wpedantic -Wconversion -Iinclude

SRCDIR = src
INCDIR = include
BINDIR = bin

all: $(BINDIR)/chain

$(BINDIR):
	mkdir -p $(BINDIR)

$(BINDIR)/chain: $(SRCDIR)/main.cpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BINDIR)/chain_tests: tests/test_chain.cpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $< -o $@

test: $(BINDIR)/chain_tests
	./$(BINDIR)/chain_tests

clean:
	rm -rf $(BINDIR) build *.paf *.svg

.PHONY: all test clean
