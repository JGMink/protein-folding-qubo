# Reference-energy investigation (2026-11-02)

This investigation uses the existing `QUBO1` files, contact matrix, penalties,
and `qubo_bench.model.validate` function. It does not alter the formulation or
the serial C++ solver. The existing benchmark metadata remains the authority
for prior reference paths and their stated certainty. The repository's
`signed20` contact table differs from the Miyazawa–Jernigan table cited by the
paper, so energies from that paper should not be imported as reference values
for these files. No published optimum has been established for these exact
exported instances and coefficient table.

## Audited instances and current records

"Lowest QUBO" means the lowest evaluated assignment known in this repository,
valid or invalid. It is an empirical upper bound on the complete QUBO minimum,
except for `toy16`. The feasible column is restricted to configurations that
pass the independent validator. A complete Hamiltonian-path enumeration proves
the 64- and 576-variable feasible minima; it does not prove their complete-QUBO
minima.

| Instance | Variables; nonzero terms | Lowest QUBO assignment known | Best validated feasible folding | Reference certainty |
|---|---:|---:|---:|---|
| `toy16` | 16; 84 | −1.00 | −1.00 | Complete-QUBO optimum by 65,536-state enumeration. |
| `fold64` | 64; 1,184 | −2.84 | −2.84 | Feasible optimum by all 28 oriented Hamiltonian paths. Complete-QUBO optimum unknown. |
| `fold576` | 576; 44,000 | −10.05 | −10.05 | Feasible optimum by all 7,220 oriented Hamiltonian paths, confirmed again here. Complete-QUBO optimum unknown. |
| `witness1024` | 1,024; 108,656 | **−17.90, invalid** (new tabu sample) | −14.25 | Valid witness only. Neither optimum is proven. |

The established feasible configurations are embedded under `reference.path` in
each benchmark file and generated from `qubo_bench.export.PRESETS`. The new
invalid −17.90 assignment and all other pilot bitstrings are in the dated JSONL
files under `experiments/`. No new feasible record was found. The incomplete
32-residue path search also saved a different validated −14.25 path.

## Bounded pilot results

The saved 500-sweep, eight-read, seed-42 baselines are
[`comparison_seed42.jsonl`](../experiments/comparison_seed42.jsonl). The dated
pilot below used seed 20261102 and the installed `dwave-samplers` 1.7.0 and
`dimod` 0.12.21. SA used eight reads of 5,000 sweeps, geometric temperature
10 to 0.1 (23.04 million attempted flips for `fold576`, 40.96 million for
`witness1024`). Tabu used eight reads with a 100 ms timeout per read and tenure
20. Sampler call times exclude QUBO parsing and BQM construction. These small
batches are diagnostic, not a reliability estimate or controlled speed ranking.

| Instance | Method | Lowest final QUBO energy | Valid final folds | Best valid energy | Sampler call time |
|---|---|---:|---:|---:|---:|
| `fold576` | SA, 5,000 sweeps | −9.72 | 1/8 | −9.72 | 0.272 s |
| `fold576` | Tabu, 100 ms/read | −8.57 | 5/8 | −8.57 | 0.852 s |
| `witness1024` | SA, 5,000 sweeps | −16.33 | 0/8 | — | 0.715 s |
| `witness1024` | Tabu, 100 ms/read | −17.90 | 0/8 | — | 0.942 s |

The 32-residue lattice search evaluated 45,054 complete paths and 41.7 million
search nodes in a 30-second limit. It found −14.25 but did not finish, so it
establishes no new bound. The same script completed the 24-residue search in
1.21 seconds and recovered all 7,220 paths and −10.05. A time-limited partial
path search explores starts in site order, so it is not an unbiased random
sample of all valid folds.

Every saved candidate is rescored against its exported QUBO and checked with
`qubo_bench.model.validate`. Each record includes the instance SHA-256, full
bitstring, path and contact energy when valid, software version, seed, settings,
and computational effort. For tabu, the time limit is clock based; the same
seed can yield different trajectories on another machine. Path enumeration is
deterministic and records a null seed.

## Reproduce

Run these commands from the repository root. Outputs can use other filenames
to retain the dated records. `paths` needs only the Python standard library;
`sample` needs the already optional `dimod` and `dwave-samplers` packages.

```sh
python -m experiments.search_references paths benchmarks/fold576.qubo --seconds 30 --output experiments/fold576_paths_repeat.jsonl
python -m experiments.search_references paths benchmarks/witness1024.qubo --seconds 30 --output experiments/witness1024_paths_repeat.jsonl
python -m experiments.search_references sample benchmarks/fold576.qubo --sampler sa --reads 8 --seed 20261102 --sweeps 5000 --output experiments/fold576_sa_repeat.jsonl
python -m experiments.search_references sample benchmarks/witness1024.qubo --sampler sa --reads 8 --seed 20261102 --sweeps 5000 --output experiments/witness1024_sa_repeat.jsonl
python -m experiments.search_references sample benchmarks/fold576.qubo --sampler tabu --reads 8 --seed 20261102 --timeout-ms 100 --tenure 20 --output experiments/fold576_tabu_repeat.jsonl
python -m experiments.search_references sample benchmarks/witness1024.qubo --sampler tabu --reads 8 --seed 20261102 --timeout-ms 100 --tenure 20 --output experiments/witness1024_tabu_repeat.jsonl
```

The baseline serial C++ solver and its saved outputs remain as documented in
[`serial_solver.md`](serial_solver.md). The dated pilot does not replace them.

## Search options and limits

| Option | Capability and cost | Recommendation |
|---|---|---|
| Longer D-Wave SA or existing C++ SA | Direct QUBO benchmark; easy to scale sweeps and reads. Neither enforces folding constraints. | Use bounded batches for raw QUBO comparisons, but measure validated-fold hit rates separately. |
| D-Wave multistart tabu | Direct QUBO heuristic with seed, tenure, initial states, reads, and per-read timeout. The pilot produced valid 24-residue folds and a lower invalid 32-residue QUBO value. | Retain as a baseline; do not infer physical improvement from its QUBO energy. |
| OpenJij SA/SQA | Maintained alternative samplers with sweep/read controls and custom schedules. Not installed in this environment. SQA requires extra tuning and runtime. | Optional independent QUBO pilot after a controlled installation; compare equal work budgets and validate every result. |
| `dimod.ExactSolver` | Enumerates binary assignments and is intended for small test models. | Appropriate for 16 variables; infeasible at 64 or more. |
| OR-Tools CP-SAT or SCIP | Time-limited constraint/MIP models can return feasible incumbents and valid objective bounds if modeled correctly. A direct 1,024-variable QUBO has 108,656 nonzero terms, so linearization and weak relaxations may be costly. Neither library is installed here. | Pilot a *separate* lattice-aware 32-residue model with explicit path constraints and exact integer contact energies. Record bound, gap, status, and independently validate the exported path. Do not treat a feasible-only status as proof. |
| Hamiltonian-path enumeration | Directly searches valid full-lattice paths. It proves the 24-residue feasible minimum rapidly; 32-residue naive enumeration exceeded 30 seconds. | Useful for reference configurations. Consider symmetry or connectivity pruning before requesting a longer complete search. |

[D-Wave SA API](https://docs.dwavequantum.com/en/latest/ocean/api_ref_samplers/generated/dwave.samplers.SimulatedAnnealingSampler.sample.html),
[D-Wave tabu API](https://docs.dwavequantum.com/en/latest/ocean/api_ref_samplers/generated/dwave.samplers.TabuSampler.sample.html),
[dimod exact-solver guidance](https://docs.dwavequantum.com/en/latest/ocean/api_ref_dimod/sampler_composites.html),
[OpenJij documentation](https://tutorial.openjij.org/intro.html), and
[OR-Tools CP-SAT time limits](https://developers.google.com/optimization/cp/cp_tasks)
document these capabilities. CP-SAT reports a best objective bound and search
status in its [solver response](https://github.com/google/or-tools/blob/stable/ortools/sat/cp_model.proto).

## What to benchmark next

`fold576` is the strongest larger **validated-feasibility** comparison: its
feasible optimum is proven, and the pilot produced valid samples at modest
cost. Report hit rate and time to reach −10.05, with the validity check applied
to each saved sample. `witness1024` is suitable now for **raw QUBO** algorithm
comparison, with −17.90 as an empirical incumbent. Its zero valid samples in
both 500- and 5,000-sweep unconstrained batches make it a poor valid-fold
success benchmark until a lattice-aware reference search or a sampler with a
measurable valid-fold hit rate is available.

The follow-up [lattice-aware CP-SAT pilot](cpsat_reference.md) used the −14.25
path as a solution hint. It recovered that witness but found no better valid
fold and left a weak lower bound. The reference investigation is complete at
the current budgets: −10.05 remains the proven 24-residue feasible minimum,
and −14.25 remains the best validated 32-residue folding without an optimum
claim. A larger QUBO MIP run or long unconstrained sampling campaign is not
justified by these pilots.
