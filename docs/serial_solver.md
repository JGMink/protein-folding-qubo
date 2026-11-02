# Standalone serial QUBO baseline

This is a conventional binary-variable Metropolis simulated annealer. The C++17
executable reads a portable QUBO file and has **no Python runtime dependency**.
Python is used only to generate protein instances, validate saved candidates,
and run optional D-Wave sampler comparisons. No external C/C++ source was copied.

## Build

On Linux or macOS, from the repository root:

```sh
mkdir -p build
c++ -std=c++17 -O2 -Wall -Wextra -o build/qubo_sa solver/qubo_sa.cpp
```

On Windows, in an MSVC Developer PowerShell or Developer Command Prompt:

```powershell
New-Item -ItemType Directory -Force build | Out-Null
cl /std:c++17 /O2 /EHsc /W4 /Febuild\qubo_sa.exe /Fobuild\qubo_sa.obj solver\qubo_sa.cpp
```

The compiled executable is ignored by Git. This repository does not include
prebuilt binaries.

## Generate QUBOs

The initial collection is generated deterministically:

```sh
python -m qubo_bench.export --all --output benchmarks
```

Export an additional 2D or 3D case without editing evaluation code:

```sh
python -m qubo_bench.export --sequence CFIKDLWV --dims 3 4 \
  --penalties 3 4 4 --model signed20 --output benchmarks/custom.qubo
python -m qubo_bench.export --sequence ENDVTMDMDPSPCLFRIHNLPRAHSFDRFGWHQFDKYHYKWKWAWAPS \
  --dims 4 4 3 --order x-fastest --output benchmarks/custom3d.qubo
```

The generator's signed20 table was extracted unchanged from the experiments
notebook. The toy table uses the tutorial's positive contact-strength convention
and negates it when constructing the QUBO. The same E1, E2, and E3 penalties are
used in both. The tests check the extracted matrix against both notebooks and
compare generated coefficients against each relevant notebook builder.

## QUBO1 input

The UTF-8 text format has one required record per line:

```text
QUBO1
n 2
offset 3
meta {"description":"example"}
q 0 0 -1
q 0 1 0.5
q 1 1 -2
```

`n` is the number of binary variables, indexed `0..n-1`. `q i j value` records
must have `0 <= i <= j < n`; omitted coefficients are zero, and duplicate
indices are errors. Each exported coefficient is written with 17 significant
digits, enough to round-trip an IEEE double. The energy is

`offset + sum_i q(i,i)*x_i + sum_(i<j) q(i,j)*x_i*x_j`.

Protein files contain one compact JSON `meta` record with the sequence,
penalties, complete contact table, lattice dimensions, numbered site
coordinates, lattice edges, site order, and variable order. Bit `i*M+u` means
residue `i` occupies site `u`, where `M` is the number of sites. The 8×4
notebook section uses column-major site numbering; its exported coordinates
and edges preserve that numbering. The C++ optimizer ignores metadata; the
Python validator uses it to check physical feasibility independently.

## Solver CLI and results

```sh
build/qubo_sa benchmarks/fold64.qubo --sweeps 1000 --trials 20 --seed 42 \
  --schedule geometric --temp-start 10 --temp-end 0.1 \
  --checkpoint 100 --target -2.84 --output build/fold64.jsonl
build/qubo_sa benchmarks/toy16.qubo --energy 0010100001000001
build/qubo_sa benchmarks/toy16.qubo --delta 0010100001000001 0
```

On Windows use `build\qubo_sa.exe`. `--help` lists all options. Defaults are
1,000 sweeps, one trial, seed 1, geometric temperature from 10 to 0.1, and a
checkpoint every 100 sweeps. One sweep attempts one flip per variable in index
order. Trials use independent PRNG streams seeded `seed + trial`. An explicit
`--initial` bitstring starts each trial at that state; otherwise each bit is
random. The `fixed` schedule uses `--temp-start` for every sweep. For linear
and geometric schedules, the first and last sweeps use the stated endpoints;
a one-sweep run uses the final temperature. Zero sweeps returns the initial
state. Temperatures must be positive.

The flip change is `(1-2*x_i) * (q(i,i) + sum_(j!=i) q(i,j)*x_j)`.
Neighbor-weighted local fields are initialized once and updated only after
accepted flips. A downhill or equal-energy flip is always accepted; an uphill
flip is accepted with probability `exp(-delta/T)`. The solver records the
initial state, best state after each accepted flip, and final state. Its JSONL
output has one object per trial with seed, settings, attempted and accepted
updates, runtime, best/final energies and bitstrings, checkpoint progression,
and optional `first_hit`. The target is a numerical QUBO energy threshold;
it is **not** a physical-feasibility test and does not stop a trial early.
Checkpoints record the initial state,
each requested sweep interval, and the final sweep, without logging every
attempt. `--checkpoint 0` leaves only initial and final records.

## Reproducible experiments

Run both samplers on the same instance with the same number of sweeps and
trials:

```sh
python -m qubo_bench.experiment benchmarks/fold64.qubo \
  --solver build/qubo_sa --neal --sweeps 500 --trials 8 --seed 42 \
  --schedule geometric --temp-start 10 --temp-end 0.1 \
  --output experiments/comparison.jsonl
```

On Windows use `--solver build/qubo_sa.exe`. The script also accepts several
instance paths and a `--neal-only` mode. It accepts `neal` if installed, or
the current `dwave.samplers` successor. The comparison maps temperature to
inverse temperature using `beta = 1/T`; both runs receive `n * sweeps`
attempts per trial. The samplers have different random streams and may have
different update details, so identical trajectories are not expected. Compare
their **final** states for like-for-like sampling; the C++ result additionally
contains its best encountered state. C++ `runtime_seconds` measures each trial
inside the executable. The Python sampler's `batch_runtime_seconds` is the
whole `sample()` call and is repeated on each read record, so these time fields
are not a controlled speed comparison.

The experiment script recomputes each saved candidate's QUBO energy and checks
physical feasibility from the metadata. Its JSONL file retains bitstrings and
paths for later independent analysis. It does not check feasibility at every
attempt. A feasible-reference success requires a valid saved incumbent and
an energy at or below the reference; it cannot be inferred from QUBO energy
alone. A run may visit a good feasible state that is not its best QUBO state,
so this measurement can undercount feasible hits.

## Initial references and measured baseline

| Instance | Variables | Reference category | Energy | Basis |
|---|---:|---|---:|---|
| `toy16` | 16 | Proven complete-QUBO global minimum | −1.00 | All 65,536 binary states enumerated with penalties (1,1,1); eight minimizing states, all valid. |
| `fold64` | 64 | Proven feasible-fold minimum | −2.84 | All 28 Hamiltonian paths on 2×4 enumerated; complete QUBO optimum unproven. |
| `fold576` | 576 | Proven feasible-fold minimum | −10.05 | All 7,220 Hamiltonian paths on 4×6 enumerated; complete QUBO optimum unproven. |
| `witness1024` | 1,024 | Known valid configuration | −14.25 | A checked self-avoiding path; no optimum claim. |

The initial collection is only a starting point. The export CLI supports other
sequences, lattice sizes, site orders, and penalty vectors. The 48-residue
4×4×3 case can be exported with the command above, without asserting the
paper's reported energy for this repository's distinct contact table.

`experiments/baseline_neal_seed42.jsonl` records eight reads per instance,
500 sweeps per read, seed 42, geometric temperatures 10 to 0.1. It was run
with `dwave-samplers` 1.7.0 and `dimod` 0.12.21. The table describes **final
samples**, not optimality guarantees:

| Instance | Best final QUBO energy | Valid final samples | Best valid contact energy |
|---|---:|---:|---:|
| `toy16` | −1.00 | 8/8 | −1.00 |
| `fold64` | −2.84 | 5/8 | −2.84 |
| `fold576` | −4.75 | 0/8 | — |
| `witness1024` | −15.50 | 0/8 | — |

The last QUBO energy is below the known valid witness energy, but its sampled
state is invalid. It is not evidence of a better physical folding.

`experiments/comparison_seed42.jsonl` contains the completed WSL 2 comparison:
eight 500-sweep reads per instance, seed 42, geometric temperatures 10 to 0.1,
using `dwave-samplers` 1.8.0 and `dimod` 0.12.23. Every saved energy was
recomputed against its exported QUBO, and physical validity was checked from
the embedded lattice metadata. These are empirical samples, not optimum
proofs. C++ final states and D-Wave final states are shown on the same basis:

| Instance | C++ best final energy | C++ valid final | D-Wave best final energy | D-Wave valid final |
|---|---:|---:|---:|---:|
| `toy16` | −1.00 | 8/8 | −1.00 | 8/8 |
| `fold64` | −2.84 | 6/8 | −2.84 | 5/8 |
| `fold576` | −4.12 | 0/8 | −4.75 | 0/8 |
| `witness1024` | −14.99 | 0/8 | −15.50 | 0/8 |

The C++ run also records the best QUBO state encountered during each read. Its
lowest incumbents were −1.00, −2.84, −4.58, and −15.10 respectively; the last
two were invalid folds. These small batches do not establish a performance
ranking. Both methods received the same nominal update budget and temperature
endpoints, while their transition order and random streams differ.

## Tests, provenance, and limitations

```sh
python -m pytest -q -p no:cacheprovider tests
QUBO_SA_EXE="$PWD/build/qubo_sa" python -m pytest -q -p no:cacheprovider tests
```

On Windows set `$env:QUBO_SA_EXE = (Resolve-Path build/qubo_sa.exe).Path` before
the second command. The tests cover serialization, exact toy references,
notebook coefficient equivalence, Python/C++ energy and flip changes,
Metropolis acceptance, schedules, seeding, zero sweeps, and larger-instance
execution. `neal` checks are optional when that package is unavailable.

On the development Windows 11 Home host, MSVC compiled the solver successfully
with `/W4`, but Smart App Control blocked the unsigned executable. The Code
Integrity Operational log recorded an enforced event 3077 under policy
`VerifiedAndReputableDesktop` (GUID `0283ac0f-fff1-49ae-ada1-8a933130cad6`).
Smart App Control's registry state was 1 (on). Freshly installed D-Wave package
DLLs were blocked too; the preinstalled `dwave.samplers` worked. Smart App
Control has no per-app allow switch, so the practical local validation path is
a Linux environment such as WSL 2, which was used for the validation below.
It can be installed without turning off Windows protection. On Windows
Home, Windows Sandbox is unavailable. See Microsoft's
[WSL install guide](https://learn.microsoft.com/en-us/windows/wsl/install) for
the administrator install and restart steps.

For local validation while keeping Smart App Control on, install WSL 2 with
Ubuntu from an administrator PowerShell (`wsl --install -d Ubuntu`) and restart
when prompted. In Ubuntu, go to this repository through `/mnt/c/...`, then:

```sh
sudo apt update
sudo apt install -y g++ python3-venv
python3 -m venv ~/protein-qubo-venv
source ~/protein-qubo-venv/bin/activate
python -m pip install numpy pytest dimod dwave-samplers
mkdir -p build
c++ -std=c++17 -O2 -Wall -Wextra -o build/qubo_sa solver/qubo_sa.cpp
QUBO_SA_EXE="$PWD/build/qubo_sa" python -m pytest -q -p no:cacheprovider tests
python -m qubo_bench.experiment benchmarks/toy16.qubo benchmarks/fold64.qubo \
  benchmarks/fold576.qubo benchmarks/witness1024.qubo \
  --solver "$PWD/build/qubo_sa" --neal --sweeps 500 --trials 8 --seed 42 \
  --output experiments/comparison_seed42.jsonl
```

The Linux executable is built and run inside WSL; the Windows host policy
stays enabled. A separate Linux machine offers the same workflow. Review the
test result and each solver's final-state distribution before treating the
comparison as validated.

On November 2, 2026, Ubuntu 26.04.1 on WSL 2 built the solver with GCC 15.2
using `-Wall -Wextra`, ran it successfully, and reported **26 passed** for the
full test suite with `QUBO_SA_EXE` set. This includes the three C++ runtime
checks previously skipped on Windows. The comparison above was run in the same
environment. Windows Smart App Control remained enabled.

The supplied 20-type matrix differs from the Miyazawa–Jernigan table cited by
the research paper, as recorded in the existing repository. This work does
not change the matrix or the protein-folding formulation. The original toy
model uses a different positive contact-strength convention. For provenance,
the solver is an original implementation informed by the conventional
Metropolis rule and by inspection of D-Wave's
[`dwave-neal`](https://github.com/dwavesystems/dwave-neal) and
[`dwave-samplers`](https://github.com/dwavesystems/dwave-samplers) projects.
Their source code was not reused or adapted, so their Apache 2.0 license does
not apply to the new C++ file. No published interaction-energy matrix was
substituted.
