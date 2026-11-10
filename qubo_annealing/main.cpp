#include "experiment.h"

// Entry point: run_experiment (experiment.cpp) does everything; its return value is the exit code (0 ok, 2 error).
// An MPI version would wrap this call in MPI_Init / MPI_Finalize.
int main(int argc, char** argv) {
    return run_experiment(argc, argv);
}
