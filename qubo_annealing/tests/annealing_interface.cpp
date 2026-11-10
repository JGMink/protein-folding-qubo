#include "annealing.h"
#include <cassert>
#include <future>
int main() {
    const Qubo problem{2, 0, {-1, 2}, {{{1, -0.5}}, {{0, -0.5}}}, {1, 2, {{1}, {0}}}};  // lattice: one residue, two adjacent sites
    AnnealingParameters parameters;
    parameters.schedule = "paper25";
    parameters.selection = "random";
    parameters.sweeps_per_temperature = 3;
    auto first = std::async(std::launch::async, [&] { return simulated_annealing(problem, parameters, 19); });
    auto second = std::async(std::launch::async, [&] { return simulated_annealing(problem, parameters, 19); });
    const auto a = first.get();
    const auto b = second.get();
    assert(a.final_state == b.final_state && a.best_state == b.best_state);
    assert(a.accepted_updates == b.accepted_updates && a.attempted_updates == 150);
    assert(a.progress.size() == 26 && a.final_energy == energy(problem, a.final_state));
    assert(a.runtime_seconds >= a.progress.back().elapsed_seconds);
    parameters.initial_state = {2, 0};
    bool rejected = false;
    try { simulated_annealing(problem, parameters, 1); }
    catch (const std::invalid_argument&) { rejected = true; }
    assert(rejected);
}
