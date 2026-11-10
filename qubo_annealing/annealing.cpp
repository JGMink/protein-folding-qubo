#include "annealing.h"
#include "fold_tracker.h"
#include <chrono>
#include <cmath>
#include <random>
#include <limits>
#include <stdexcept>
using Clock = std::chrono::steady_clock;

// Draw only for uphill moves, preserving the seeded proposal stream.
static bool accept_metropolis_move(double delta, double temperature, std::mt19937_64& rng) {
    if (delta <= 0) return true;
    return std::generate_canonical<double, 53>(rng) < std::exp(-delta / temperature);
}

// Give this sweep's temperature, including the one-sweep endpoint convention.
static double temperature_for_sweep(const AnnealingParameters& parameters, std::uint64_t sweep) {
    if (parameters.schedule == "paper25")
        return std::pow(1.05, -double(1 + sweep / parameters.sweeps_per_temperature));
    if (parameters.schedule == "fixed") return parameters.start;
    double fraction = parameters.sweeps <= 1 ? 1.0 : double(sweep) / double(parameters.sweeps - 1);
    if (parameters.schedule == "linear") return parameters.start + (parameters.end - parameters.start) * fraction;
    return parameters.start * std::pow(parameters.end / parameters.start, fraction);
}

// Run one Metropolis trial against a read-only QUBO. Each trial owns its random
// generator, configuration, cached fields, and incumbent; a sweep visits each
// variable once in index order by default, or draws with replacement when
// random selection is requested. Return best and final states separately.
AnnealingResult simulated_annealing(const Qubo& problem, const AnnealingParameters& requested_parameters, std::uint64_t trial_seed) {
    AnnealingParameters parameters = requested_parameters;
    if (problem.n < 1 || problem.diagonal.size() != static_cast<std::size_t>(problem.n) || problem.neighbors.size() != static_cast<std::size_t>(problem.n))
        throw std::invalid_argument("invalid QUBO dimensions");
    if (parameters.start <= 0 || parameters.end <= 0 || !std::isfinite(parameters.start) || !std::isfinite(parameters.end) ||
        (parameters.schedule != "linear" && parameters.schedule != "geometric" && parameters.schedule != "fixed" && parameters.schedule != "paper25") ||
        (parameters.selection != "sequential" && parameters.selection != "random") ||
        (parameters.schedule == "paper25" && parameters.sweeps_per_temperature == 0))
        throw std::invalid_argument("invalid annealing parameters");
    if (parameters.schedule == "paper25") {
        if (parameters.sweeps_per_temperature > std::numeric_limits<std::uint64_t>::max() / 25)
            throw std::invalid_argument("too many sweeps per temperature");
        parameters.sweeps = 25 * parameters.sweeps_per_temperature;
        parameters.checkpoint = parameters.sweeps_per_temperature;
    }
    if (parameters.has_target && !std::isfinite(parameters.target)) throw std::invalid_argument("invalid target");
    // A supplied fold layout must describe this QUBO, or the validity checks would index out of bounds.
    if (parameters.fold && (parameters.fold->residues < 1 || parameters.fold->sites < 1 ||
        static_cast<long long>(parameters.fold->residues) * parameters.fold->sites != problem.n ||
        parameters.fold->next_to.size() != static_cast<std::size_t>(parameters.fold->sites)))
        throw std::invalid_argument("fold layout does not match the QUBO");
    if (!parameters.initial_state.empty()) {
        if (parameters.initial_state.size() != static_cast<std::size_t>(problem.n)) throw std::invalid_argument("initial state has wrong size");
        for (unsigned char bit : parameters.initial_state)
            if (bit > 1) throw std::invalid_argument("initial state must be binary");
    }
    AnnealingResult result{};
    const auto started = Clock::now();
    const auto elapsed = [&]() {
        return std::chrono::duration<double>(Clock::now() - started).count();
    };
    result.seed = trial_seed;
    std::mt19937_64 rng(trial_seed);
    std::uniform_int_distribution<int> random_bit(0, problem.n - 1);

    // An explicit initial bitstring is reused for each trial; otherwise draw one
    // from this trial's seeded generator.
    std::vector<unsigned char> state(problem.n);
    if (!parameters.initial_state.empty()) state = parameters.initial_state;
    else for (auto& bit : state) bit = static_cast<unsigned char>(rng() & 1);

    // fields[bit_index] is the selected bit's coefficient given the other current bits. Building
    // all fields once lets each proposed flip use a local delta instead of a
    // full QUBO energy calculation; accepted flips update only affected fields.
    auto fields = local_fields(problem, state);
    double current_energy = energy(problem, state), best_energy = current_energy;
    auto best_state = state;
    std::uint64_t attempted_updates = 0, accepted_updates = 0;

    // The initial state may already be the run's best or meet the numerical
    // energy target, so record sweep/update zero before any proposals.
    result.first_best = {0, 0, elapsed()};
    result.target_reached = parameters.has_target && best_energy <= parameters.target + 1e-9;
    if (result.target_reached) result.first_hit = {0, 0, elapsed()};
    result.progress.push_back({0, 0, best_energy, elapsed(), current_energy, {}});

    // Optional fold measurements (see fold_tracker.h); they never touch the RNG.
    FoldTracker tracker;
    if (parameters.fold) {
        tracker = fold_tracker_start(*parameters.fold, state);
        if (fold_beats_record(tracker, current_energy) && fold_is_valid_now(tracker, state))  // the start may already be valid
            fold_record(tracker, problem, state, {0, 0, elapsed()});
        fold_fill_validity(tracker, state, result.progress.front());  // the sweep-0 sample pushed above
    }

    for (std::uint64_t sweep = 0; sweep < parameters.sweeps; ++sweep) {
        double current_temperature = temperature_for_sweep(parameters, sweep);
        for (int attempt = 0; attempt < problem.n; ++attempt) {
            const int bit_index = parameters.selection == "random" ? random_bit(rng) : attempt;
            double delta = flip_delta(state, fields, bit_index);
            ++attempted_updates;

            // Metropolis accepts downhill/equal flips. For an uphill flip,
            // reject when a uniform draw is at least exp(-delta / T).
            if (!accept_metropolis_move(delta, current_temperature, rng))
                continue;
            ++accepted_updates;

            // change is +1 for 0->1 and -1 for 1->0. Update neighboring fields
            // before changing the bit, then add the same delta to current energy.
            int change = 1 - 2 * int(state[bit_index]);
            update_fields_after_flip(problem, fields, bit_index, change);
            state[bit_index] ^= 1;
            current_energy += delta;

            // Every strict improvement replaces the incumbent. The last such
            // event is the first encounter of the run's eventual best energy.
            if (current_energy < best_energy - 1e-9) {
                best_energy = current_energy; best_state = state;
                result.first_best = {sweep + 1, attempted_updates, elapsed()};
            }

            // A target is only a QUBO energy threshold. Record its first hit
            // independently of the best state, and continue the anneal.
            if (!result.target_reached && parameters.has_target && current_energy <= parameters.target + 1e-9) {
                result.target_reached = true;
                result.first_hit = {sweep + 1, attempted_updates, elapsed()};
            }

            // Fold measurements: keep the tracker's counts current, then log a new record-low valid fold.
            if (parameters.fold) {
                fold_update_counts(tracker, bit_index, change);
                if (fold_beats_record(tracker, current_energy) && fold_is_valid_now(tracker, state))
                    fold_record(tracker, problem, state, {sweep + 1, attempted_updates, elapsed()});
            }
        }

        if (parameters.fold) fold_end_of_sweep(tracker, state, current_energy);  // per-sweep sample for the checkpoint stats

        // Store only requested checkpoints and the final sweep, not each flip.
        if ((parameters.checkpoint && (sweep + 1) % parameters.checkpoint == 0) || sweep + 1 == parameters.sweeps) {
            result.progress.push_back({sweep + 1, attempted_updates, best_energy, elapsed(), current_energy, {}});
            if (parameters.schedule == "paper25") {
                auto& point = result.progress.back();
                point.current_energy = energy(problem, state);
                for (auto bit : state) point.bits += char('0' + bit);
            }
            if (parameters.fold)
                fold_fill_checkpoint(tracker, result.progress.back(), problem, state, attempted_updates, accepted_updates);
        }
    }

    result.attempted_updates = attempted_updates;
    result.accepted_updates = accepted_updates;
    result.best_state = best_state;
    result.final_state = state;

    // Incremental energy can accumulate rounding error; evaluate saved states
    // once in full before reporting their energies.
    result.best_energy = energy(problem, best_state);
    result.final_energy = energy(problem, state);
    if (parameters.fold) {
        result.best_is_valid = is_valid_fold(*parameters.fold, best_state);
        result.final_is_valid = is_valid_fold(*parameters.fold, state);
        result.valid_history = tracker.history;
        result.best_valid_state = tracker.best_state;
    }
    result.runtime_seconds = elapsed();
    return result;
}
