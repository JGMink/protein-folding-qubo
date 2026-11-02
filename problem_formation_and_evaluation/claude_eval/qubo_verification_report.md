# Historical verification note

The November 2025 review formerly in this file concluded that a few hand-picked
energy equalities verified the entire implementation. That conclusion was too broad.
The checks also assumed the number of residues equaled the number of lattice sites,
used a case-sensitive import path incorrectly, and did not require each expected
component or feasibility status to match.

The current checks are in `tests/test_qubo_correctness.py`. They include exhaustive
binary-state evaluation for the four-residue toy model, independent constraint
calculations, sign checks for each experiments-notebook QUBO builder, and
Hamiltonian-path benchmarks for the supplied contact matrix. The notebook's
20-type matrix does not match the Miyazawa–Jernigan table cited by
Irbäck et al.; its benchmarks apply to the supplied matrix and 2D lattices.
