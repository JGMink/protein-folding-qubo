# Serial QUBO simulated annealing

A C++17 solver with no Python runtime dependency. Start here; see the
[developer walkthrough](../docs/serial_readiness.md), [MJ-48 audit](../docs/irback2025_audit.md),
and [curated papers](../docs/serial_references.md) for details.

| File | Responsibility |
|---|---|
| `serial.cpp` | Short executable entry point |
| `experiment.h/.cpp` | CLI validation, diagnostics, trials, JSONL output |
| `annealing.h/.cpp` | Independent parameters, result, and Metropolis annealing |
| `qubo.h/.cpp` | Sparse QUBO loading, energy, local fields, flip updates |
| `Makefile` | C++17 optimized build |
| `tests/test_solver.py` | Arithmetic, reproducibility, witnesses, protocol checks |
| `run/run_experiment.sh` | One protected parameter-driven run |
| `run/prepare_mj48.sh` | Selected or all six MJ inputs with published penalties |
| `run/` | Also retains historical smoke and audit scripts |
| `benchmarks/` | Fixed fold8/fold24/fold32 inputs |
| `analysis/summarize.py` | Summarize experiment logs |

From the repository root under Linux or WSL (C++17, make, Bash, Python 3):

```sh
make -C hpc_serial
python3 -m unittest discover -s hpc_serial/tests -v
# Quick workflow check: 128 proposals, no research run.
bash hpc_serial/run/run_experiment.sh fold8 smoke-check --sweeps 2 --seed 7 --checkpoint 1
# Select a small benchmark and ordinary solver parameters.
bash hpc_serial/run/run_experiment.sh fold24 small-example --schedule geometric --temp-start 570 --temp-end 0.00885 --sweeps 10 --trials 2 --seed 42
# Prepare a selected MJ-48 input, with published penalties (1.5, 2, 2).
bash hpc_serial/run/prepare_mj48.sh mj5
# Cheap protocol example: 25 stages x 1 sweep, not a paper-length experiment.
bash hpc_serial/run/run_experiment.sh mj5 mj-example --schedule paper25 --selection random --sweeps-per-temperature 1 --seed 42
# Summarize exactly the selected logs; CSV goes to standard output.
python3 hpc_serial/analysis/summarize.py hpc_serial/results/smoke-check/trials.jsonl hpc_serial/results/small-example/trials.jsonl > hpc_serial/results/example-summary.csv
```

These are usage examples, not a final experiment matrix or meaningful performance
measurements. Use a fresh run identifier each time. The runner refuses an existing
directory, even for a failed run; absence of `status.txt` containing `complete`
means the run did not finish successfully. Results live in `results/RUN_ID/`:
`trials.jsonl`, one-argument-per-line `arguments.txt`, UTC start time, input/binary
SHA-256 checksums, and a Makefile snapshot. Checksums identify exact artifacts;
the Makefile snapshot does not capture compiler versions or overridden build flags.
Preserve the corresponding source and executable when archiving an experiment.

Set `RESULTS_DIR=/path/to/results` before the runner to choose another parent.
The first argument accepts `fold8`, `fold24`, `fold32`, `mj1` through `mj6`, or an
existing QUBO path (relative to your current directory). `MJ_DIR` overrides the MJ
input directory for both scripts. Run `prepare_mj48.sh` without arguments to
prepare all six, or pass selected names. It refuses existing files. Python is
needed for generation; it is not needed by the solver or single-run script.
The runner never rewrites input coefficients. MJ aliases should point to inputs
created by the preparation script; arbitrary supplied files retain their own model.

All remaining runner arguments pass directly to `serial`; use `hpc_serial/serial --help` for the complete list. `--output` is managed by the runner. `--initial BITS`
sets an initial binary configuration, `--target ENERGY` requests target tracking,
and `--trials N` uses independent RNG seeds `seed, seed+1, ...`. JSONL stores the
zero-based trial index and actual seed. Keep seeds within the uint64 range,
including increments for repetitions. A direct solver command remains supported:

```sh
hpc_serial/serial hpc_serial/benchmarks/fold8.qubo --sweeps 2 --seed 7
```

Direct `--output` overwrites its destination; use the runner for protected outputs.
For parameter lists, ordinary shell loops suffice; give every setting a distinct
run ID. Runs execute serially:

```sh
for seed in 7 8; do
    bash hpc_serial/run/run_experiment.sh fold8 "example-seed-$seed" --sweeps 2 --seed "$seed"
done
```

The summarizer accepts explicit file lists or shell-expanded globs such as
`hpc_serial/results/example-seed-*/trials.jsonl`. It groups by input path,
benchmark, schedule, selection, sweeps/stage, temperatures, total sweeps,
checkpoint interval, initial configuration, and target. It reports trial counts,
best and final energy statistics, mean runtime, and numeric target-hit rates.
Supply each log once; repeated/overlapping seeds are counted as supplied, not deduplicated.
Target-hit rates are not physical-fold success rates. Different paths to identical
inputs remain separate groups; changed files at the same path are not detected by
the summarizer, so consult the run checksums.

The fixed 8/24/32-residue examples are development benchmarks, distinct from
paper-compatible MJ-48 (48 residues, 2,304 bits, full 4x4x3 grid). The historical
`run_serial.sh` and `run_irback2025_pilot.sh` retain their fixed lists and overwrite
named outputs; they are not the standard entry points for new experiments.
`run_smoke.sh` similarly regenerates its dedicated smoke fixture.

Best energy/configuration means the incumbent retained during the trajectory;
final means the state after all updates. They can differ. Reported state energies
are recomputed in full, including the QUBO offset. A numeric target is an energy
threshold, not a feasibility or optimality certificate; hitting it does not stop a run.

One sweep attempts N flips. Sequential selection visits each bit once; random
selection draws N times with replacement. Linear/geometric temperatures change
per sweep; fixed uses the starting temperature. `paper25` uses 25 plateaus with
T=1.05^(-stage), stages 1-25, overriding total sweeps and checkpoint spacing.
Trial seeds are seed + trial index. Checkpoints record progress; paper stage
checkpoints additionally save current bits/energy and beta. `first_best` and
`first_hit` record the first relevant sweep/update counts. Runtime covers algorithm
initialization, updates, checkpoints, and final energy recomputation, ending before
JSON serialization; loading and output are excluded.

Preliminary verification checked all 210 MJ matrix entries against the source,
12 native assignments, 24 invalid assignments, and eight saved MJ pilot histories.
All six native witness energies agree with published values. The pilot reached
MJ5's native energy but does not reproduce published SA success rates: two seeds
and at most 10,000 sweeps/stage versus 100 trials and up to 80,000 in the paper.
Witnesses establish attainability; global optimality relies on published enumeration.
Feasibility was checked for saved states, not every transient state. Exact author
RNG/update-order equivalence remains unknown. See the audit for results and provenance.

Workflow verification: seven solver/workflow tests pass, including exclusive run
creation, seed/trial recording, and separate summary groups for update modes.
A generated MJ5 input was checked for penalties (1.5, 2, 2); its inexpensive
25-sweep run recorded 57,600 proposals and the initial plus 25 stage checkpoints.
These checks validate execution and logging, not research performance.
