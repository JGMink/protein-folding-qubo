#ifndef HPC_SERIAL_ANNEALING_H
#define HPC_SERIAL_ANNEALING_H
#include "qubo.h"
#include <cstdint>

struct AnnealingParameters {
    std::string schedule = "geometric";
    std::string selection = "sequential";
    std::vector<unsigned char> initial_state;
    std::uint64_t sweeps = 1000, sweeps_per_temperature = 1000, checkpoint = 100;
    double start = 10, end = 0.1, target = 0;
    bool has_target = false;
};
// Where the final best energy or a requested target was first reached.
struct Event {
    std::uint64_t sweep = 0, attempted_updates = 0;
    double elapsed_seconds = 0;
};

// A compact progress sample taken at the start and at selected sweep boundaries.
struct Checkpoint {
    std::uint64_t sweep, attempted_updates;
    double best_energy, elapsed_seconds;
    double current_energy = 0;
    std::string bits;
};

// Everything produced by one anneal; serialization happens after the run.
struct AnnealingResult {
    std::uint64_t seed, attempted_updates, accepted_updates;
    std::vector<unsigned char> best_state, final_state;
    double best_energy, final_energy;
    Event first_best, first_hit;
    bool target_reached;
    std::vector<Checkpoint> progress;
    double runtime_seconds = 0;
};

// Each call owns its RNG and working state; the problem is read-only.
AnnealingResult simulated_annealing(const Qubo& problem, const AnnealingParameters& parameters, std::uint64_t seed);
#endif
