CXX = c++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

.PHONY: all clean check-interface
all: qubo_anneal

qubo_anneal: main.cpp experiment.cpp experiment.h annealing.cpp annealing.h fold_tracker.cpp fold_tracker.h qubo.cpp qubo.h
	$(CXX) $(CXXFLAGS) -o $@ main.cpp experiment.cpp annealing.cpp fold_tracker.cpp qubo.cpp

clean:
	rm -f qubo_anneal

# Link the reusable algorithm without the CLI, and exercise independent callers.
check-interface:
	mkdir -p ../build/readiness
	$(CXX) $(CXXFLAGS) -pthread -I. tests/annealing_interface.cpp annealing.cpp fold_tracker.cpp qubo.cpp -o ../build/readiness/interface_test
	../build/readiness/interface_test
