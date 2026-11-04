# Independent MJ-48 audit and bounded serial reproduction

Audit date: 2026-11-03 (America/New_York).

**Recommendation: GO for beginning MPI parallel-tempering development on the
verified numerical model.** The source-derived witnesses, interaction matrix,
constraint polynomial, and serial arithmetic pass the checks below. **Published
SA success rates have not been reproduced:** this is an eight-trial pilot, not
a statistical replication. Keep this serial protocol as the baseline before
making performance claims. No MPI/PT was implemented, and nothing was committed
or pushed. Existing uncommitted work was preserved.

## Scope and initial findings

Inspected `docs/mj85_benchmarks.md`, `qubo_bench/{mj85_benchmarks.py,
mj85_table_vi.json,model.py,export.py,io.py,experiment.py}`, the MJ85 tests,
`solver/qubo_sa.cpp`, and `qubo_annealing/{main.cpp,qubo.cpp,qubo.h,README.md,
run/run_serial.sh,run/run_smoke.sh,tests/test_solver.py}`.

Before changes or pilots, the audit plan was: compare primary-source equations
and all matrix entries; recover independent reference structures if possible;
make only opt-in protocol changes; then run MJ5/MJ1 with short published budgets.
The initial discrepancies were doubled penalties, continuously changing sweep
temperatures, sequential variable order, unrelated fixed-input run scripts, and
the need to compare final states separately from incumbents.

The existing tests already supported export round trips, symmetry, checksums,
geometry, sparse-QUBO arithmetic, incremental deltas, and seeded reproducibility.
Duplicated sequences/reference values, a checksum of a local transcription, and
two scorers reading that same transcription did **not** independently establish
source accuracy or reference optimality. `test_matrix_full_audit_and_default_unchanged`
has a stronger-sounding name than its internal consistency assertions justify.
The independent evidence added here addresses those gaps without relabeling the
old tests as external validation.

## Paper versus implementation

Primary source: [Irbäck, Knuthson, Mohanty (2025), Sections II A–F,
III A and Appendices A–B](https://juser.fz-juelich.de/record/1047045/files/8n7p-7lh2.pdf).

| Item | Paper | Original repository | Audit comparison |
|---|---|---|---|
| Instances | Six MJ-48 sequences | MJ presets present; standalone run script uses different 8/24/32-residue inputs | MJ5 and MJ1 presets |
| Lattice | Open, fully occupied 4×4×3 | Matching MJ preset | Unchanged; 2,304 variables |
| Interactions | 1985 MJ Table VI, physical term written −C | Signed Table VI energies | epsilon = −C; once per nonbonded contact |
| Penalty weights | (1.5,2,2) | (3,4,4) in MJ export defaults | Explicit (1.5,2,2); global defaults preserved |
| Initialization | Random bits | Independent seeded random bits unless `--initial` | Random bits, no warm start |
| Move/acceptance | Single-bit Metropolis | Matching cached-field deltas | Unchanged |
| Temperatures | 25 plateaus, beta_i=1.05^i, i=1..25 | Linear/geometric temperature changes every sweep, or fixed | Opt-in `paper25`, T_i=1/beta_i |
| Bit selection | Not explicitly specified; sweep is an average of one attempt per variable | Sequential | Opt-in uniform random with replacement; explicit assumption |
| Work | 1,000–80,000 sweeps **per temperature** | `--sweeps` is total | 1,000 and 10,000 per stage; 25,000/250,000 total |
| Trials | 100 per setting | Standalone script defaults to five on other inputs | Two seeds (42,43) per setting |
| Stopping | Specified run length | Continues after numeric target hit | Full budget, no early stop or target supplied |
| Comparison statistic | Final energy and returned-solution success | Saves incumbent and final separately | Both verified/reported; only final matches are comparable successes |

Paper temperature endpoints are approximately 0.952381 and 0.295303.
Both C++ implementations originally use sequential selection and sweep-wise
temperature schedules. Only `qubo_annealing/qubo_anneal.cpp` receives the opt-in changes;
`solver/qubo_sa.cpp` and existing scripts retain their behavior. Choosing random
selection is plausible given the paper's sweep wording, but exact author-code
equivalence cannot be certified without their update-order details/RNG/code.

## Physical model and matrix evidence

All six sequences and energies were compared character by character with
2025 Table I and independently with the image of [Faísca and Plaxco (2006),
Table 2](https://cdn.ncbi.nlm.nih.gov/pmc/blobs/6eb0/2242573/e3987467b879/1608tbl2.jpg).
No differences were found. Topology 1 corresponds to A (MJ1–3), topology 2 to
B (MJ4–6).

The original [Miyazawa–Jernigan (1985) scan, printed p.545, Table VI](https://www.sanzo.org/~miyazawa/profile/.p/publications/85_Macromolecules_18_534-552/fulltext.pdf)
was downloaded afresh, PDF page 12 rendered and rotated 90 degrees, and all
210 lower-triangle coefficients visually transcribed again into
`experiments/irback2025_audit/matrix_scan_transcription.json`. **210/210 match**
the repository. This is a fresh source comparison by the audit, not a second
human annotator or an author-provided machine-readable table. Source download
hashes are retained for provenance, not presented as correctness evidence.

Residue order is `CMFILVWYAGTSQNEDHRKP`, with Trp before Tyr. Values are signed
contact-preference energies at printed precision 0.01, in dimensionless RT
units. The caption identifies the transformation e_ij+e_rr−e_ir−e_jr; Table VI
is already transformed, so applying that transformation again is wrong.
Anchors include C–C=−1.06, M–W=−0.67, E–K=−0.97, and P–P=+0.26.
No rescaling, offset, sign reversal, or replacement with Table V/1996 MJ values
is needed. All six externally recovered native scores agreeing to the printed
hundredth strongly support this interpretation of the 2025 authors' −C notation.
Their exact machine-readable matrix has not been obtained.

Geometrically there are 48 sites and 104 undirected face-neighbor edges:
36+36+32. Each valid fold has 47 covalent edges and hence **57 nonbonded
contacts**. Contacts exclude adjacent sequence indices and are counted once.
Full occupancy follows from 48 singly placed residues without overlaps on
48 sites; an extra vacancy constraint is unnecessary.

The generator's polynomial agrees with the intended per-contact interpretation
of Eqs. (3)–(6):

- E1=sum_i(sum_n b_i,n−1)^2 expands to offset 48 lambda1, diagonal −lambda1,
  and pair coefficient 2 lambda1 within each residue row.
- E2=one half of the ordered i≠j sum equals the implemented i<j overlap pairs,
  with coefficient lambda2 (no extra factor of two).
- E3 penalizes consecutive residues at Euclidean distance strictly greater than
  one. On this integer grid that means distinct, nonadjacent sites. Same-site
  consecutive placements are excluded from E3 and penalized by E2, as in the paper.
- The contact polynomial uses i<j−1 and both orientations of each undirected
  lattice edge. A particular occupied nonbonded contact contributes one epsilon.

`experiments/irback2025_audit/verify.py` calculates equations directly from
coordinate sets and the fresh scan transcription, without the builder's edge
list, table loader, or validator. It checks all twelve native QUBO assignments
(six sequences × two penalty settings), plus 24 random/empty/dense invalid
assignments against the expanded QUBO. These invalid checks cover the polynomial
away from the feasible manifold; they are not a proof of its global minimum.

Doubling penalties changes an arbitrary state's energy by
1.5 E1+2 E2+2 E3 relative to the paper. For a lone residue-count defect that is
1.5 extra; for a lone overlap/connectivity defect, 2 extra. It leaves valid-fold
energies unchanged but changes barriers and Metropolis probabilities.
It is not equivalent to scaling beta because the interaction term is unchanged.
The default offsets are 144 versus the paper-weight offset of 72. Include the
offset when comparing full QUBO energy; omitting it shifts every score by 72 here.
No global guarantee that invalid states cannot beat a fold follows merely from
positive penalties, source agreement, or these finite tests.

## Genuine source-derived ground-state witnesses

The 2025 data-availability statement says data are available from the authors
on reasonable request. No 2025 coordinate file or supplement was found in the
paper, linked publication records, or focused author-resource search. The public
[author repository](https://github.com/sandipan-mohanty/DWaveHPLatticeProteins)
is for the older 2D HP work, not these MJ-48 data. No authors were contacted.

However, [2006 Figure 1](https://cdn.ncbi.nlm.nih.gov/pmc/blobs/6eb0/2242573/4f683f5e00c9/1608fig1.jpg)
contains genuine native **contact maps** for both target topologies. In the
720×599 source image, the two plots span approximately x=54..341 and
x=428..714, y=258..544, on axes 0..50. Sampling integer residue positions
1..48 (black grayscale <100), then converting to zero-based indices, gives
114 ordered marks per topology: exactly 57 symmetric, noncovalent pairs.
There were no missing/asymmetric marks or energy-driven edits to the maps.
Adding the 47 chain bonds gives 104-edge graphs isomorphic to a 4×4×3 grid.
NetworkX graph isomorphism reconstructs coordinates **without using sequence
energies**. The stored coordinates are this audit's reconstruction from a
published contact map, not author-supplied coordinates. Rotations/reflections
are immaterial.

| Benchmark | Published energy | Independent contact sum | Full QUBO, either penalty setting | E1,E2,E3 |
|---|---:|---:|---:|---|
| MJ1 | −25.85 | −25.85 | −25.85 | 0,0,0 |
| MJ2 | −25.92 | −25.92 | −25.92 | 0,0,0 |
| MJ3 | −26.09 | −26.09 | −26.09 | 0,0,0 |
| MJ4 | −25.87 | −25.87 | −25.87 | 0,0,0 |
| MJ5 | −26.15 | −26.15 | −26.15 | 0,0,0 |
| MJ6 | −26.24 | −26.24 | −26.24 | 0,0,0 |

These independently establish **attainability** and the exact energies of
the reported native structures. Global optimality still rests on the authors'
published enumeration, not a new enumeration in this repository. A stored
reference number alone was never used as a witness.

## Actual bounded results

Random initialization, random selection with replacement, paper penalties,
25 beta plateaus, no tuning between runs. Every best/final/stage score was
recomputed by the separate equation scorer. Runtime is the solver's trial
timer (initial fields and checkpoint/output overhead included; QUBO loading
excluded), under WSL Ubuntu, GCC 15.2.0, C++17 `-O2`, Ryzen 7 7800X3D.
The two instances were launched as independent processes with overlapping
execution; these are not controlled single-process runtime benchmarks.

| Instance | Sweeps/stage | Seed | Best full QUBO | Best valid? / physical | Final full QUBO | Final valid? | Seconds |
|---|---:|---:|---:|---|---:|---|---:|
| MJ5 | 1,000 | 42 | −26.15 | yes / −26.15 | −26.15 | yes | 1.32 |
| MJ5 | 1,000 | 43 | −17.64 | no / undefined | −16.17 | no | 1.35 |
| MJ1 | 1,000 | 42 | −16.86 | no / undefined | −8.42 | no | 1.38 |
| MJ1 | 1,000 | 43 | −16.94 | no / undefined | −15.92 | no | 1.41 |
| MJ5 | 10,000 | 42 | −26.15 | yes / −26.15 | −23.29 | no | 13.32 |
| MJ5 | 10,000 | 43 | −26.15 | yes / −26.15 | −26.15 | yes | 14.11 |
| MJ1 | 10,000 | 42 | −20.25 | yes / −20.25 | −18.38 | no | 14.52 |
| MJ1 | 10,000 | 43 | −19.68 | no / undefined | −19.67 | no | 13.80 |

Per trial, attempted updates are 57,600,000 or 576,000,000 respectively;
total 2,534,400,000. Final ground-state matches are MJ5 1/2 at each budget,
MJ1 0/2. Best-visited ground-state matches are MJ5 1/2 then 2/2, MJ1 0/2.
For MJ5 seed 42 at 10,000 sweeps/stage, the first optimum incumbent occurs
in stage 19 (sweep 181,820), and stage 20 is a valid optimum, yet the final
state has one missing/extra-row defect (E1=1) and energy −23.29. This is why
incumbent success cannot substitute for final-state success.

Example stage-end **current** full-QUBO energies (not monotonically decreasing;
Metropolis accepts uphill moves):

| Instance, budget, seed | Stage 1 | Stage 5 | Stage 10 | Stage 15 | Stage 20 | Stage 25 |
|---|---:|---:|---:|---:|---:|---:|
| MJ5, 1,000, 42 | 47.74 | 51.92 | 29.96 | 13.00 | −6.93 | −26.15 |
| MJ5, 10,000, 42 | 58.84 | 48.38 | 26.60 | 17.38 | −26.15 | −23.29 |
| MJ1, 10,000, 42 | 63.35 | 41.17 | 22.38 | 12.07 | −5.72 | −18.38 |

All 25 stage states, decomposed interactions/penalties, accepted-update counts,
and incumbent histories are in `verified_summary.json` and the four raw JSONL
logs. Invalid-state interaction energies are available there but deliberately
not labeled physical-fold energies. Absence of a valid saved state does not
prove that no transient valid state was visited: feasibility is checked at
saved incumbents, final states, and stage boundaries, not every attempted flip.

The pilot exhibits energy/constraint improvement and source-energy attainment
for MJ5, qualitatively compatible with the paper. Its failure to solve MJ1
at these budgets does not reproduce or refute the longest-budget result.
The paper's strongest SA claims concern 80,000 sweeps/stage and 100 trials;
the pilot uses at most one eighth that per-trial work and two trials. We did
not digitize short-budget success probabilities from the paper's plots and
cannot make a calibrated agreement claim. No parameter tuning was attempted.

Local mean sweep cost is approximately 0.055 ms versus the paper's reported
approximately 1.1 ms; cached fields, hardware, implementation, environment and
timing scope differ. This is not evidence of matching dynamics or a measured
speedup against the authors' code. A 100-trial 80,000/stage setting would require
460.8 billion proposals; extrapolating locally gives roughly three hours per
sequence, versus roughly 58 hours using the paper's reported long-run timing.
No such run was launched.

## Validation, limitations, and next step

- `python -m pytest -q -p no:cacheprovider --basetemp build/audit/pytest-mj tests/test_mj85_benchmarks.py`: **17 passed**. The first attempt had temp-directory permission errors; rerunning with a workspace temp directory resolved them.
- `python3 -m unittest discover -s qubo_annealing/tests -v` in WSL: **6 passed**, including existing delta/energy/reproducibility checks, staged update counts/endpoints, checkpoint energies, and random-with-replacement behavior.
- `python experiments/irback2025_audit/verify.py`: **210 matrix entries, 12 native assignments, 24 invalid assignments, and all eight saved trial histories passed**. The direct scorer uses the independent scan transcription.
- `python3 experiments/irback2025_audit/check_cpp.py` in WSL: native full-size MJ1/MJ5 C++ energies and twelve flip deltas agreed with full QUBO recomputation.
- The final standalone build passed `-Wall -Wextra` with no warnings. The earlier pilot binary differed only in explicit initialization of newly added checkpoint fields; both binary hashes are recorded.

There is no evidence from these checks of a matrix, geometry, contact-count,
penalty-polynomial, or incremental-energy bug. Original defaults and schedules
were unsuitable for a faithful paper comparison; the opt-in experiment fixes
that discrepancy. Remaining uncertainty concerns author update ordering/RNG,
their exact matrix/code, finite-budget trapping, and statistical performance.
The data do not isolate the cause of remaining SA failures; computational budget
and stochastic trapping are plausible, not established diagnoses. No claim that
larger penalties explain a measured failure is made, since no controlled
penalty-ablation experiment was run.

Begin MPI PT against this verified model and retained serial baseline. Before
claiming paper-equivalent SA success, obtain the authors' bit-selection details
and run a separately costed larger seed set/longer budget with **final-state
feasibility** as the success criterion. Optional future improvements include
tracking the best feasible state during annealing and incremental penalty counts;
neither is required to validate the present energy function, and neither was
silently added to improve pilot success.

## Reproduce and inspect

From the repository root, `bash qubo_annealing/run/run_irback2025_pilot.sh` builds
and runs the bounded pilot, regenerating its named logs and verified summary.
The standard-library verification script also runs independently on saved logs.
Source downloads/renderings are in ignored `build/audit`; source URLs and hashes
are retained below that experiment directory in `source_hashes.json`.
`digitize.py` requires Pillow and the downloaded Figure 1 at
`build/audit/native_maps.jpg`; `reconstruct.py` additionally needs NetworkX,
reads the extracted map, and reconstructs coordinates without energy input.
Its fallback package search includes the repository's existing `build/pylibs`.
The frozen contact/coordinate fixture means these extra dependencies are not
required to verify all reported energies. Witness contact indices are zero-based.
The original global benchmark metadata was left unchanged; the witness fixture
is explicitly separate and source-provenanced by this report.

Path note: the serial implementation directory is now named `qubo_annealing/`. The historical
audit manifest retains its original paths and hashes as recorded at audit time.
