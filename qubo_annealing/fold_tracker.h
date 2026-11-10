#ifndef QUBO_FOLD_TRACKER_H
#define QUBO_FOLD_TRACKER_H
#include "annealing.h"

// Fold bookkeeping for one annealing chain, kept apart from the algorithm. It watches the chain's
// state flip by flip, logs each new record-low valid fold, and gathers per-checkpoint statistics.
// It never touches the random generator, so it cannot change the trajectory.
struct FoldTracker {
    const FoldLayout* layout = nullptr;     // lattice description (not owned)
    std::vector<int> ones_in_residue;       // per residue: how many of its bits are 1
    int placed = 0;                         // how many residues have exactly one 1-bit
    std::vector<ValidRecord> history;       // each new lowest valid energy, oldest first
    std::vector<unsigned char> best_state;  // lowest-energy valid fold seen (empty if none)
    // Statistics since the previous checkpoint:
    std::uint64_t block_sweeps = 0, block_valid_sweeps = 0;
    double block_energy_sum = 0;            // sum of the end-of-sweep energies
    std::uint64_t attempts_at_last_checkpoint = 0, accepts_at_last_checkpoint = 0;
};

// Start tracking from `state`. `layout` must outlive the tracker.
FoldTracker fold_tracker_start(const FoldLayout& layout, const std::vector<unsigned char>& state);

// Call after every accepted flip of bit `bit_index` (change = +1 for 0->1, -1 for 1->0).
void fold_note_flip(FoldTracker& tracker, int bit_index, int change);

// Is `state` a valid fold right now? Cheap test first; the O(n) scan only runs if every residue is placed.
bool fold_is_valid_now(const FoldTracker& tracker, const std::vector<unsigned char>& state);

// True if `state` (running energy `energy_now`) is a valid fold lower than every valid fold recorded so far.
bool fold_is_new_record(const FoldTracker& tracker, const std::vector<unsigned char>& state, double energy_now);

// Log `state` as the new record, stamped with `when`. The energy is recomputed in full, so it has no rounding drift.
void fold_record(FoldTracker& tracker, const Qubo& problem, const std::vector<unsigned char>& state, const Event& when);

// Call at the end of every sweep: adds this sweep's ending energy and validity to the current block.
void fold_end_of_sweep(FoldTracker& tracker, const std::vector<unsigned char>& state, double energy_now);

// Fill `point.valid` and the best-valid fields from the tracker (this alone is enough for the sweep-0 sample).
void fold_fill_validity(const FoldTracker& tracker, const std::vector<unsigned char>& state, Checkpoint& point);

// Fill all fold fields of a checkpoint, including the block statistics, then start a fresh block.
void fold_fill_checkpoint(FoldTracker& tracker, Checkpoint& point, const Qubo& problem,
                          const std::vector<unsigned char>& state,
                          std::uint64_t attempted_updates, std::uint64_t accepted_updates);

#endif
