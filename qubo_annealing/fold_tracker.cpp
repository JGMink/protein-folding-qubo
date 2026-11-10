#include "fold_tracker.h"

// A valid fold has exactly one 1-bit in each residue's block of `sites` bits, so we count those
// per residue. That count lets most validity checks finish without scanning the whole state.
FoldTracker fold_tracker_start(const FoldLayout& layout, const std::vector<unsigned char>& state) {
    FoldTracker tracker;
    tracker.layout = &layout;
    tracker.ones_in_residue.assign(layout.residues, 0);
    for (std::size_t bit = 0; bit < state.size(); ++bit)  // bit / sites is the residue this bit belongs to
        tracker.ones_in_residue[bit / layout.sites] += state[bit];
    for (int count : tracker.ones_in_residue) tracker.placed += (count == 1);
    return tracker;
}

// Only the flipped bit's residue changes, so updating `placed` is O(1).
void fold_update_counts(FoldTracker& tracker, int bit_index, int change) {
    int& count = tracker.ones_in_residue[bit_index / tracker.layout->sites];  // the flipped bit's residue
    tracker.placed -= (count == 1);  // leaving or entering "exactly one" changes `placed`
    count += change;
    tracker.placed += (count == 1);
}

// The count test is free; the O(n) scan is only reached when every residue is placed exactly once.
bool fold_is_valid_now(const FoldTracker& tracker, const std::vector<unsigned char>& state) {
    return tracker.placed == tracker.layout->residues  // the counts cover only rule 1: one site per residue
           && is_valid_fold(*tracker.layout, state);   // the scan also checks shared sites and chain adjacency
}

// "Lower" means by more than 1e-9, so floating-point noise never counts as a new record.
bool fold_beats_record(const FoldTracker& tracker, double energy_now) {
    return tracker.history.empty() || energy_now < tracker.history.back().energy - 1e-9;
}

void fold_record(FoldTracker& tracker, const Qubo& problem, const std::vector<unsigned char>& state, const Event& when) {
    tracker.history.push_back({energy(problem, state), when});  // exact energy: no rounding drift
    tracker.best_state = state;                                 // copies every bit
}

void fold_end_of_sweep(FoldTracker& tracker, const std::vector<unsigned char>& state, double energy_now) {
    ++tracker.block_sweeps;
    tracker.block_energy_sum += energy_now;
    if (fold_is_valid_now(tracker, state)) ++tracker.block_valid_sweeps;
}

void fold_fill_validity(const FoldTracker& tracker, const std::vector<unsigned char>& state, Checkpoint& point) {
    point.valid = fold_is_valid_now(tracker, state);
    point.has_best_valid = !tracker.history.empty();
    if (point.has_best_valid) point.best_valid_energy = tracker.history.back().energy;
}

// Block statistics cover the sweeps since the previous checkpoint (or since the start).
void fold_fill_checkpoint(FoldTracker& tracker, Checkpoint& point, const Qubo& problem,
                          const std::vector<unsigned char>& state,
                          std::uint64_t attempted_updates, std::uint64_t accepted_updates) {
    point.current_energy = energy(problem, state);  // full recompute: no rounding drift
    fold_fill_validity(tracker, state, point);
    point.has_block = true;
    point.acceptance_rate = double(accepted_updates - tracker.accepts_at_last_checkpoint) /
                            double(attempted_updates - tracker.attempts_at_last_checkpoint);
    point.mean_energy = tracker.block_energy_sum / double(tracker.block_sweeps);
    point.valid_fraction = double(tracker.block_valid_sweeps) / double(tracker.block_sweeps);
    tracker.block_sweeps = tracker.block_valid_sweeps = 0;  // start a fresh block
    tracker.block_energy_sum = 0;
    tracker.attempts_at_last_checkpoint = attempted_updates;
    tracker.accepts_at_last_checkpoint = accepted_updates;
}
