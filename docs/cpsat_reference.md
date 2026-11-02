# Lattice-aware CP-SAT reference search

This is a separate reference search for the exported 32-residue, 8×4 instance.
It uses the same signed 20-type interaction coefficients and lattice metadata
as `benchmarks/witness1024.qubo`. It does not change the QUBO, its penalty
weights, the serial SA solver, or the established reference witness.

## Model and fit

The 32 residues fill all 32 sites. `position[i]` gives residue `i`'s site and
`occupant[u]` gives the residue at site `u`. CP-SAT's inverse constraint makes
these arrays permutations, so every residue has one site and sites cannot
overlap. Consecutive `position` variables must form one of the directed lattice
edges, enforcing chain connectivity. For each of the 52 physical lattice edges,
an element lookup returns the interaction energy of its endpoint residues;
chain bonds contribute zero. The 52 costs are summed and minimized.

All contact values are scaled to exact integer hundredths before entering
CP-SAT. Every reported incumbent, including intermediate solutions, is
converted to the repository's residue-major QUBO bitstring and checked by
`qubo_bench.model.validate` and `qubo_bench.model.energy`. The saved objective
bound is a **lower bound on the feasible-fold energy**, not on the minimum of
the unconstrained QUBO. The existing −14.25 path is supplied as a search hint.

This is a reasonable compact model: 64 site-index variables, 31 adjacency
tables, and 52 contact lookups. Building it took about 0.01–0.02 seconds in
the pilots. The difficult part is proving a good bound: the edge-cost lookups
do not give CP-SAT a strong relaxation of the global path energy.

## Bounded pilots (2026-11-02)

OR-Tools 9.15.6755 was installed in an isolated temporary directory for these
runs. Both used seed 20261102 and the original `witness1024.qubo` file. Solver
time excludes Python startup and model construction. Four-worker runs can vary
slightly between machines despite a fixed seed.

| Instance | Workers | Limit | Status | Best validated fold | Best objective bound | Validated incumbents |
|---|---:|---:|---|---:|---:|---:|
| 32 residues, 8×4 | 1 | 20 s | FEASIBLE | −14.25 | −57.20 | 1 |
| 32 residues, 8×4 | 4 | 30 s | FEASIBLE | −14.25 | −52.80 | 2 |
| 8 residues, 2×4 (control) | 1 | 5 s | OPTIMAL | −2.84 | −2.84 | 1 |

The 32-residue paths returned by CP-SAT are independently valid, but neither
run improved the existing witness. Their wide gaps provide no useful evidence
that −14.25 is optimal. `FEASIBLE` means only that a fold was found before the
limit; `OPTIMAL` on the control instance means CP-SAT closed the gap for that
separate model. The JSON records preserve paths, full bitstrings, energies,
timing, status, bound, seed, worker count, software version, and the SHA-256 of
the benchmark file.

## Reproduce

Use a separate virtual environment so the optional OR-Tools dependencies do
not alter the repository's existing sampler environment:

```powershell
python -m venv .venv-cpsat
.\.venv-cpsat\Scripts\python.exe -m pip install ortools==9.15.6755
.\.venv-cpsat\Scripts\python.exe -m experiments.cpsat_reference benchmarks/witness1024.qubo --seconds 20 --workers 1 --seed 20261102 --output experiments/witness1024_cpsat_repeat.json
```

On Linux or macOS, use `.venv-cpsat/bin/python` instead. The saved pilot files
are `experiments/witness1024_cpsat_pilot_20261102.json`,
`experiments/witness1024_cpsat_pilot4_20261102.json`, and
`experiments/fold64_cpsat_check_20261102.json`.

## Recommendation

The reference investigation is complete at the current budgets. The 20- and
30-second CP-SAT pilots made no feasible-energy progress and left very loose
bounds. The existing time-limited Hamiltonian-path enumerator in
`experiments/search_references.py` remains available for reproducibility; it
reproduced the proven 24-residue minimum but did not finish the 32-residue
search in 30 seconds. A certified 32-residue minimum would require a
substantially stronger exact model or a faster complete path enumeration, so
no longer search is planned from these findings.

OR-Tools documents [time limits](https://developers.google.com/optimization/cp/cp_tasks)
and the [solver status and best objective bound](https://github.com/google/or-tools/blob/stable/ortools/sat/cp_model.proto).
