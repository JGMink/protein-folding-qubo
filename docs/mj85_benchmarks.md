# MJ85 Table VI and the six MJ-48 benchmarks

**Audit update:** [Independent audit and bounded SA reproduction](irback2025_audit.md)
supersedes the original integration limitations below. Figure 1 of Faísca and
Plaxco (2006) supplies native contact maps for both topologies. The audit
digitized them, reconstructed lattice coordinates, and independently recovered
all six reference energies, including zero-penalty complete QUBO assignments.
The general benchmark metadata still leaves `ground_state_coordinates=None`;
the source-derived witnesses live in `experiments/irback2025_audit/native_folds.json`.
Existing checksum and shared-table tests are internal consistency gates, not
independent source verification. A fresh scan transcription and independent
equation scorer are retained with the audit.

The six 48-residue instances are defined in `qubo_bench/mj85_benchmarks.py`.
Their sequences, topology labels, and reported minima were checked against
[Table I of Irbäck, Knuthson, and Mohanty (2025)](https://juser.fz-juelich.de/record/1047045/files/8n7p-7lh2.pdf)
and the original sequence image in
[Table 2 of Faísca and Plaxco (2006)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2242573/).
The `reported_exhaustively_established` status describes the published finite-grid
enumeration. This repository has not reproduced those minima, and indexed ground
state coordinates/contact maps were unavailable for this integration.

## Interaction resource

`qubo_bench/mj85_table_vi.json` is the separately selectable, version 1
transcription of the signed contact-preference energies from
[Miyazawa and Jernigan (1985), Table VI, page 545](https://www.sanzo.org/~miyazawa/profile/.p/publications/85_Macromolecules_18_534-552/fulltext.pdf).
It stores the 210 printed values in integer hundredths of RT, in the printed
lower-triangular arrangement. The loader expands them to a symmetric 20 × 20
matrix with 400 populated cells. The SHA-256 of the JSON bytes is recorded in
`qubo_bench/mj85_table_vi.sha256` and checked at load time.

The printed residue order is **CYS MET PHE ILE LEU VAL TRP TYR ALA GLY THR
SER GLN ASN GLU ASP HIS ARG LYS PRO**, or `CMFILVWYAGTSQNEDHRKP`.
The implementation handoff reversed Trp and Tyr; the scan places **Trp before
Tyr**. The 210 entries were transcribed visually from the original PDF scan
after rotating page 12, then checked visually a second time row by row. The
source anchors are Cys–Cys = −1.06 and Pro–Pro = +0.26. Selected checks using
the printed Table V relation, `epsilon_ij = e_ij + e_rr - e_ir - e_jr`, give
Cys–Cys = −1.06, Pro–Pro = +0.26, and Met–Trp = −0.67. Published intermediate
terms are rounded, so some other reconstructed entries differ by 0.01.

These are **signed energies** and are summed once per noncovalent contact. The
2025 paper writes its physical term as `-sum C`; when using the original signed
Table VI directly, `epsilon = -C`. No sign reversal, offset, or normalization
is applied to `MJ85_TableVI`. The existing `signed20` default and `toy`
resources are unchanged. The authors' exact machine-readable matrix has not
been obtained, so the equality of every numerical coefficient to their code
remains unverified.

## Geometry and scoring

Every benchmark uses a 4 × 4 × 3 open simple-cubic grid, all 48 sites occupied,
and 47 face-adjacent covalent bonds. The grid has
`(4−1)×4×3 + 4×(4−1)×3 + 4×4×(3−1) = 104` undirected edges. Therefore every
valid Hamiltonian path has `104 − 47 = 57` distinct nonbonded contacts.

`evaluate_coordinates` validates bounds, occupancy, continuity, and contact
count without using QUBO metadata or its edge list. It sums the physical
energy in integer hundredths. `build` and `validate` use the existing sparse
QUBO formulation, where the offset is `lambda_1 × 48`. For a valid assignment,
all three constraint penalties vanish and the offset-corrected QUBO energy
equals the independent physical score. Invalid assignments cannot match a
published reference: experiment records require physical validity and exact
agreement at two-decimal precision; a score below the reference is flagged for
investigation.

## Export and verify

```powershell
python -m qubo_bench.export --preset irback2025_mj5 --output benchmarks/irback2025_mj5.qubo
python -m qubo_bench.export --preset irback2025_mj1 --output benchmarks/irback2025_mj1.qubo
python -m pytest -q -p no:cacheprovider tests/test_mj85_benchmarks.py
```

The same `--preset` interface accepts `irback2025_mj2` through
`irback2025_mj6`. The default penalties remain `(3, 4, 4)` and can be supplied
explicitly with `--penalties`. `--all` still exports only the original four
repository instances. The MJ5 pilot export has 2,304 variables and 425,248
nonzero upper-triangular terms (about 8.5 MB on disk). The independent test
fold scores −1.19 RT with 57 contacts; it is a validation fixture, not a
ground state. No comparative optimizer sweeps were run for this integration.
