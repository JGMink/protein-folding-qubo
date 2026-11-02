# 🧬 Protein Folding via QUBO Formulation

An **educational implementation** and explanation of the QUBO-based protein folding approach from:

> **"Folding lattice proteins confined on minimal grids using a quantum-inspired encoding"**  
> *Anders Irbäck, Lucas Knuthson, and Sandipan Mohanty (2025)*  
> Physical Review E 112, 045302  
> arXiv:2510.01890 · [arXiv abstract](https://arxiv.org/abs/2510.01890) · [DOI: 10.1103/8n7p-7lh2](https://doi.org/10.1103/8n7p-7lh2)

This project provides a **pedagogical walkthrough** of their methodology using **toy examples** to make the concepts accessible and verifiable.

---

## 🎯 Purpose of This Project

This repository is designed to:

1. **Explain** the QUBO formulation from Irbäck et al. (2025) with detailed derivations
2. **Demonstrate** the approach using small, tractable examples (4-6 residues)
3. **Validate** that the implementation correctly encodes the problem
4. **Provide** tested, reusable code for educational purposes

**This is NOT a production protein folding solver** - it's an educational tool to understand how to formulate biological optimization problems for quantum computing.

---

## 📄 Original Paper Summary

The paper by Irbäck, Knuthson, and Mohanty presents a hybrid quantum-classical approach that:
- Formulates lattice protein folding as a QUBO problem
- Uses a paper-inspired contact-energy objective; the supplied 20-type coefficient table does not match the Miyazawa-Jernigan Table VI cited by the paper
- Encodes physical constraints (connectivity, self-avoidance) as penalty terms
- Can be solved using quantum annealers, QAOA, or classical methods

**Our contribution**: We provide a clean, well-documented implementation with toy examples that make the mathematics accessible to students and researchers new to quantum optimization.

---

## 🌱 What This Implementation Does

This project systematically demonstrates the QUBO formulation through:

1. **Problem Definition**: Simplified 2D lattice model with 3 amino acid types (H, P, C)
2. **Mathematical Formulation**: Complete derivation of all energy terms from the paper
3. **QUBO Construction**: Step-by-step matrix building with visualization
4. **Validation**: Dual calculation methods verify correctness
5. **Toy Examples**: Small chains (4-6 residues) that can be verified by hand

The tutorial uses an illustrative three-type contact matrix. The experiments notebook uses a separate signed 20-type matrix. Both encode the same lattice constraints, but their contact coefficients have opposite sign conventions before they enter the energy function.

---

## 📊 Project Structure
```
protein-folding-qubo/
├── problem_formation_and_evaluation/
│   ├── energy_calc/
│   │   └── calc_mods.py              # Direct energy calculation
│   ├── QUBO_construction/
│   │   ├── qubo_generation.py        # QUBO matrix building
│   │   └── construction_test.py      # Validation tests
│   └── claude_eval/
│       └── verify_qubo.py            # Comprehensive verification
├── notebooks/
│   └── protein_folding_qubo.ipynb    # Main tutorial notebook
├── papers/
│   └── Simulated_and_Hybrid_Annealing_Algorithms...pdf  # Our own annealing writeup
└── README.md
```

---

## 🔍 Energy Formulation (Following Irbäck et al., 2025)

The total energy combines four terms:

$$E_{\text{total}} = E_{\text{MJ}} + \lambda_1 E_1 + \lambda_2 E_2 + \lambda_3 E_3$$

Where:
- **E_MJ**: Nonbonded nearest-neighbor contact energy (using the supplied interaction coefficients)
- **E₁**: One position per residue constraint
- **E₂**: Self-avoidance constraint (no overlaps)
- **E₃**: Chain connectivity constraint

### Binary Encoding

Each configuration uses binary variables:
- $b_{i,n} = 1$ if residue $i$ occupies position $n$
- $b_{i,n} = 0$ otherwise

For N residues and M positions: **N × M binary variables**

### QUBO Matrix Representation

All terms combine into a single matrix **Q**:

$$E(x) = x^T Q x + \text{offset}$$

This formulation is compatible with quantum annealers (D-Wave), QAOA, and classical solvers.

---

## ⚙️ Getting Started

### Requirements
```bash
pip install -r requirements.txt
```

### Running the Tutorial Notebook

#### Option 1: Local Jupyter
```bash
git clone <your-repo-url>
cd protein-folding-qubo
jupyter notebook notebooks/protein_folding_qubo.ipynb
```

#### Option 2: Google Colab
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/19EnI3qWr388mPwWXmZY4o0jzo4_UETY3?usp=sharing)

---

## 🧪 What's in the Tutorial Notebook

The notebook provides a **complete walkthrough** of the QUBO formulation:

### Section 1: Problem Definition (3-4 min)
- Introduction to the paper's approach
- Lattice protein folding basics
- Amino acid types and interactions

### Section 2: Mathematical Formulation (5-6 min)
- Binary variable encoding
- **Detailed derivation** of all four energy terms
- Step-by-step expansion of constraint penalties

### Section 3: QUBO Construction (8-10 min)
- Matrix building for each energy component
- Visualization of QUBO structure
- How constraints become penalty terms

### Section 4: Toy Example - HPCH Chain (4-5 min)
- **4 residues** on a 2×2 lattice (16 variables)
- Small enough to verify calculations by hand
- Complete matrix walkthrough

### Section 5: Validation (5-6 min)
- Valid vs. invalid configurations
- Verification: QUBO method = Direct calculation
- Constraint violation detection

### Section 6: Slightly Larger Example - HHCHPC (3-4 min)
- **6 residues** demonstrating scalability
- Still tractable for educational purposes
- Complexity analysis

### Section 7: Discussion (3-4 min)
- Scaling behavior
- Connection to quantum/classical solvers
- Extensions and future work

**Total time: ~30-35 minutes**

---

## 🧩 Validated Test Cases

All examples verified with dual calculation methods:

| Test | Chain | Variables | Status | Valid Config? |
|------|-------|-----------|--------|---------------|
| 1 | HPCH | 16 | ✓ Pass | Yes (E=−1) |
| 2 | HCHP | 16 | ✓ Pass | No (violations) |
| 3 | HHCHPC | 36 | ✓ Pass | No (E₃ violation) |

These examples check agreement between QUBO and direct energy calculations. The regression suite also exhaustively checks all 65,536 binary assignments of the four-residue toy model against an independent calculation.

---

## 🔧 Running Verification Tests

To run the regression checks:

```bash
python -m pytest -q -p no:cacheprovider tests
```

To run the standalone example checks:

```bash
cd problem_formation_and_evaluation/QUBO_construction
python construction_test.py
```

For comprehensive verification:
```bash
cd problem_formation_and_evaluation/claude_eval
python verify_qubo.py
```

Passing example checks establish agreement for those configurations; they do not establish that a chosen penalty vector makes every QUBO minimum physically valid.

---

## 📊 Visualization Features

The notebook includes:
- **QUBO matrix heatmaps**: See the structure of Q
- **Lattice diagrams**: 2D grid with residue placements
- **Energy breakdowns**: Contribution from each term
- **Color-coded violations**: Visual constraint checking

---

## 🎯 Potential Use Cases (Following the Paper)

This formulation can be used with:

1. **Quantum Annealers** (D-Wave): Direct QUBO implementation
2. **QAOA**: Convert to Ising Hamiltonian for gate-based quantum computers
3. **Simulated Annealing**: Classical probabilistic optimization
4. **Hybrid Solvers**: Quantum-classical approaches (as in the original paper)

---

## 🧠 Educational Goals

### What You'll Learn:
- How to formulate constrained optimization as QUBO
- Converting biological constraints to penalty terms
- Practical QUBO matrix construction
- Verification techniques for quantum-compatible formulations

### Why Small Examples?
- **Verifiable**: Can check calculations by hand
- **Transparent**: Every matrix entry has clear meaning
- **Debuggable**: Easy to trace errors
- **Educational**: Focus on concepts, not scale

**For production protein folding**, see the original paper's full implementation and benchmarks.

---

## 📚 Related Resources

### Original Paper
- Irbäck, A., Knuthson, L., & Mohanty, S. (2025). "Folding lattice proteins confined on minimal grids using a quantum-inspired encoding." *Physical Review E* 112, 045302.
- arXiv:2510.01890 · [arXiv abstract](https://arxiv.org/abs/2510.01890) · [DOI: 10.1103/8n7p-7lh2](https://doi.org/10.1103/8n7p-7lh2)

> **Note on the paper PDF:** we don't redistribute the publisher's typeset PDF here (it's copyrighted by the American Physical Society) — use the links above to read it.

### Background Reading
- **QUBO formulations**: Andrew Lucas (2014) - "Ising formulations of many NP problems"
- **HP model**: Ken Dill's work on simplified protein models
- **Miyazawa-Jernigan**: The cited folding paper identifies Table VI of the 1985 contact-potential paper as its coefficient source; the supplied notebook coefficients differ from that table
- **QAOA**: Farhi et al. (2014) - Original quantum approximate optimization paper

---

## 🔄 Differences from Original Paper

| Aspect | Original Paper | This Implementation |
|--------|----------------|---------------------|
| Scale | Production-scale chains | Toy examples (4-6 residues) |
| Purpose | Research results | Educational tool |
| Focus | Performance benchmarks | Mathematical clarity |
| Code | Research codebase | Pedagogical, documented |
| Validation | Quantum hardware results | Dual calculation verification |

**This is complementary to the paper**: We provide the educational foundation to understand their methodology.

---

## 🏡 Future Extensions

Potential enhancements (maintaining educational focus):

- [ ] Integration with D-Wave Ocean SDK (live quantum annealing)
- [ ] QAOA implementation with Qiskit
- [ ] Step-by-step solver comparison (SA vs exact vs quantum)
- [ ] Interactive visualization (plotly/dash)
- [ ] Additional small examples with detailed walkthroughs
- [ ] Jupyter widgets for parameter exploration

---

## 📧 Contact & Attribution

**Implementation by:** Jonah Minkoff  
**Email:** jonah@planetminkoff.com

**Based on the work of:**  
Irbäck, A., Knuthson, L., & Mohanty, S. (2025) - "Folding lattice proteins confined on minimal grids using a quantum-inspired encoding"  
*Physical Review E* 112, 045302 · arXiv:2510.01890

---

## 📜 Citation

If you use this educational implementation, please cite both:

**This implementation:**
```bibtex
@misc{minkoff2025qubo,
  author = {Minkoff, Jonah},
  title = {Protein Folding via QUBO Formulation: An Educational Implementation},
  year = {2025},
  url = {<your-repo-url>}
}
```

**Original paper:**
```bibtex
@article{irback2025folding,
  title={Folding lattice proteins confined on minimal grids using a quantum-inspired encoding},
  author={Irb{\"a}ck, Anders and Knuthson, Lucas and Mohanty, Sandipan},
  journal={Physical Review E},
  volume={112},
  pages={045302},
  year={2025},
  doi={10.1103/8n7p-7lh2}
}
```

---

## ⚖️ License

This implementation is licensed under the MIT License. See [LICENSE.txt](LICENSE.txt) for details.

**Please cite the original paper** when using this methodology:
- Irbäck, A., Knuthson, L., & Mohanty, S. (2025). "Folding lattice proteins confined on minimal grids using a quantum-inspired encoding." *Physical Review E* 112, 045302. arXiv:2510.01890
```

---

**Thanks for exploring quantum optimization for biology!** 🚀  

If this helps your understanding, please:
- ⭐ Star the repo
- 📖 Read the original paper
- 🔗 Share with others learning quantum computing or computational biology
