CXX = c++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

.PHONY: all clean check-interface
all: serial

serial: serial.cpp experiment.cpp experiment.h annealing.cpp annealing.h qubo.cpp qubo.h
	$(CXX) $(CXXFLAGS) -o $@ serial.cpp experiment.cpp annealing.cpp qubo.cpp

clean:
	rm -f serial

# Link the reusable algorithm without the CLI, and exercise independent callers.
check-interface:
	mkdir -p ../build/readiness
	$(CXX) $(CXXFLAGS) -pthread -I. tests/annealing_interface.cpp annealing.cpp qubo.cpp -o ../build/readiness/interface_test
	../build/readiness/interface_test
