#ifndef QUBO_ANNEALING_H
#define QUBO_ANNEALING_H
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

// A valid fold that set a new record low among valid folds, and when it happened.
struct ValidRecord {
    double energy = 0;
    Event when;
};

// A compact progress sample taken at the start and at selected sweep boundaries.
struct Checkpoint {
    std::uint64_t sweep, attempted_updates;
    double best_energy, elapsed_seconds;
    double current_energy = 0;
    std::string bits;
    // The fields below are the fold measurements.
    bool valid = false;            // is the current state a valid fold?
    bool has_best_valid = false;   // has any valid fold been seen yet?
    double best_valid_energy = 0;  // lowest energy among valid folds so far (if has_best_valid)
    bool has_block = false;        // false for the sweep-0 sample: no sweeps behind it yet
    double acceptance_rate = 0;    // accepted / attempted flips since the previous checkpoint
    double mean_energy = 0;        // mean of the end-of-sweep energies since the previous checkpoint
    double valid_fraction = 0;     // share of those end-of-sweep states that were valid folds
};

// Everything produced by one anneal; serialization happens after the run.
struct AnnealingResult {
    std::uint64_t seed, attempted_updates, accepted_updates;
    std::vector<unsigned char> best_state, final_state;
    double best_energy, final_energy;
    Event first_best, first_hit;
    bool target_reached;
    // Fold measurements, using the lattice stored in the QUBO.
    bool best_is_valid = false, final_is_valid = false;  // are best_state / final_state valid folds?
    // Each new lowest valid energy, oldest first: front() = first valid fold of any energy, back() = lowest and when
    // first reached. Energies strictly decrease, so the first entry at or below a threshold is its first passage.
    std::vector<ValidRecord> valid_history;
    std::vector<unsigned char> best_valid_state;         // lowest-energy valid fold seen (empty if none)
    std::vector<Checkpoint> progress;
    double runtime_seconds = 0;
};

// Each call owns its RNG and working state; the problem is read-only.
AnnealingResult simulated_annealing(const Qubo& problem, const AnnealingParameters& parameters, std::uint64_t seed);
#endif
