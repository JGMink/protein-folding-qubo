# References

**Lucas** — "Ising formulations of many NP problems." *Frontiers in Physics*, 2014.
Shows how a huge class of combinatorial problems (graph coloring, TSP, Sudoku, etc.) can be rewritten as Ising/QUBO problems, making them solvable by annealing hardware. The theoretical backbone for why QUBO is a universal formulation.

---

**Glover et al.** — "Quantum Bridge Analytics I: QUBO Formulations." *4OR*, 2022.
A practical guide to constructing QUBO matrices — how to encode objectives, constraints, and penalties. Essential reference for the mechanics of building the Q matrix used throughout this project.

---

**Kirkpatrick et al.** — "Optimization by Simulated Annealing." *Science*, 1983.
The original paper introducing simulated annealing: probabilistic hill-climbing inspired by physical annealing, where temperature controls the acceptance of worse solutions to escape local minima. The classical algorithm this project implements and extends.

---

**Farhi et al.** — "A Quantum Approximate Optimization Algorithm." *arXiv:1411.4028*, 2014.
Introduces QAOA, a gate-based quantum circuit approach to combinatorial optimization. Provides the quantum counterpart to simulated annealing and motivates why QUBO formulations matter for near-term quantum devices.

---

**D-Wave Systems** — "Leap Hybrid BQM Solver." *docs.dwavesys.com*, 2023.
Documentation for the hybrid classical-quantum solver used in experiments. The Leap hybrid solver handles problem sizes too large for the QPU alone by combining quantum annealing with classical heuristics.

---

**Preskill** — "Quantum Computing in the NISQ Era." *Quantum*, 2018.
Coins the term NISQ (Noisy Intermediate-Scale Quantum) and contextualizes the realistic capabilities and limitations of today's quantum hardware. Frames why hybrid classical-quantum methods are the practical near-term approach.

---

**Lewis & Thompson** — "Graph coloring & Sudoku." *Intl. J. Comput. Sci.*, 2022.
Demonstrates graph coloring as a concrete QUBO use case, with Sudoku as an intuitive instance. Used here to ground the abstract QUBO framework in a familiar constraint-satisfaction problem.

---

**Slepoy et al.** — "QA across 5,000 QUBO instances." *Quantum Mach. Intell.*, 2024.
Large-scale empirical benchmarking of quantum annealing across thousands of QUBO problem instances. Provides evidence for where quantum annealing outperforms, matches, or falls short of classical solvers in practice.

---

**Irbäck, Knuthson & Mohanty** — "Folding lattice proteins on minimal grids using quantum-inspired encoding." *Phys. Rev. E*, 2025.
The primary foundational paper this research builds from. A prime example of a grid-based QUBO approach to protein folding — their encoding strategy and formulation structure were the direct inspiration for the QUBO construction used in this project.
