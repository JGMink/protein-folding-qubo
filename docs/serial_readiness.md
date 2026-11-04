# Serial solver developer walkthrough

The maintained modular solver is `hpc_serial`; `solver/qubo_sa.cpp` and
`docs/serial_solver.md` describe the older standalone implementation.

`main()` calls `run_experiment()`. The experiment layer parses CLI arguments,
loads a QUBO once, handles energy/delta diagnostics, constructs algorithm parameters,
and calls `simulated_annealing()` for each seed. It serializes returned results.
No algorithm code depends on command-line parsing, paths, JSON, or output streams.

```cpp
#include "annealing.h"
Qubo problem = load_qubo("hpc_serial/benchmarks/fold8.qubo");
AnnealingParameters parameters;
parameters.sweeps = 1000;
AnnealingResult result = simulated_annealing(problem, parameters, 42);
```

Link only `annealing.cpp` and `qubo.cpp` to use the algorithm independently.
Parameters use familiar strings for schedule/selection and a byte vector for an
optional initial state (empty requests seeded initialization). `paper25` normalizes
total sweeps and checkpoint spacing inside the algorithm as well as in the CLI.
The existing `start` and `end` members are temperatures, not inverse temperatures.

The QUBO stores diagonal coefficients and sparse neighbor lists. Full energy is
`offset + sum(diagonal[i]*bit[i]) + sum(i<j, coefficient[i,j]*bit[i]*bit[j])`.
A local field is `diagonal[i] + sum(j, coefficient[i,j]*bit[j])`; flipping bit i
changes energy by `(1 - 2*bit[i])*field[i]`. An accepted flip updates neighboring
fields before toggling the bit, then adds this delta. Full evaluation occurs at
initialization, saved paper checkpoints, and for returned best/final states.
Sparse storage, neighbor traversal, arithmetic order, and RNG draws are preserved.

The algorithm reads the shared QUBO and owns its RNG, state, fields, incumbent,
and progress records. Mathematical utilities have no mutable globals or static
work buffers. Callers must provide a structurally valid QUBO (valid neighbor
indices, symmetric off-diagonal entries) and correctly sized binary states to
low-level utilities; the loader establishes these invariants. Callers must not
mutate a shared QUBO while another call reads it. No parallel execution is implemented.

The ZIP dated 2026-11-02 is a repository/source snapshot, containing existing
README/reference/developer documents and a background annealing PDF. It was
inspected as context, not imported as a second documentation tree. The current
post-audit files supersede that snapshot for numerical results.

The baseline executable/source and regression comparison script are retained in
ignored `build/readiness/`. Timing fields are recursively removed before comparing
entire JSON records, including stage states, seeds, counts, events, and targets.

## Refactor verification (2026-11-03)

- C++17 `-O2 -Wall -Wextra` build: no warnings.
- Existing solver correctness suite: 6 passed; smoke run passed.
- MJ benchmark Python tests: 17 passed (Windows Python; WSL lacks pytest).
- Independent equation audit: 210 matrix entries, 12 native assignments,
  24 invalid assignments, and eight retained pilot histories passed.
- Full-size MJ1/MJ5 native C++ energy and twelve flip deltas passed.
- Independent algorithm interface test links without experiment/main and runs
  two concurrent calls against the same problem; invalid initial state rejected.
  Run `make -C hpc_serial check-interface` on Linux/WSL.
- Regression matrix: 24 small-input cases spanning all four schedules, both
  selection modes, three starting seeds, two trials, checkpoints and targets;
  four MJ1/MJ5 settings spanning 1,000/10,000 sweeps per stage and two trials.
  Entire deterministic JSON records agree with the preserved baseline, recursively
  excluding runtime/elapsed fields. This checks 56 trial records in total.

The existing pilot logs remain historical audit evidence. Newly timed comparisons
run in memory and do not replace them. Runtime measurements are not controlled
performance benchmarks. Independent calls have no shared mutable implementation
state; this test is a practical reentrancy check, not a race-detector proof.
