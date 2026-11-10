CXX = c++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

# Source lists, defined once. ALGORITHM is what runs without the command-line layer.
ALGORITHM = annealing.cpp fold_tracker.cpp qubo.cpp
SOURCES = main.cpp experiment.cpp $(ALGORITHM)

.PHONY: all clean check-interface
all: qubo_anneal

qubo_anneal: $(SOURCES) $(wildcard *.h)
	$(CXX) $(CXXFLAGS) -o $@ $(SOURCES)

clean:
	rm -f qubo_anneal

# Link the reusable algorithm without the CLI, and exercise independent callers.
check-interface:
	mkdir -p ../build/readiness
	$(CXX) $(CXXFLAGS) -pthread -I. tests/annealing_interface.cpp $(ALGORITHM) -o ../build/readiness/interface_test
	../build/readiness/interface_test
