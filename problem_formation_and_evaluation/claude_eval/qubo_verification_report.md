# QUBO Protein Folding Implementation: Code Review & Verification Report

**Date:** November 14, 2025  
**Reviewer:** Claude  
**Status:** ✓ ALL TESTS PASSED - IMPLEMENTATION VERIFIED AS CORRECT

---

## Executive Summary

I have thoroughly reviewed your QUBO formulation for protein folding on a lattice. After comprehensive testing across multiple scenarios, I can confirm that **your implementation is mathematically correct and properly validated**. All energy components match between the QUBO matrix representation and direct calculation methods.

---

## Files Reviewed

1. **qubo_generation.py** - QUBO matrix construction
2. **calc_mods.py** - Energy calculation and validation
3. **protein_folding_qubo.py** - Complete notebook implementation

---

## Verification Results

### Test Suite Results

| Test Case | Description | Status |
|-----------|-------------|--------|
| Test 1 | Valid configuration (HPCH) | ✓ PASS |
| Test 2 | Invalid configuration with E3 violation (HHCHPC) | ✓ PASS |
| Test 3 | Multiple constraint violations (HCHP) | ✓ PASS |

### Component-Level Verification

All four energy components were verified independently:

| Component | Formula | QUBO Implementation | Direct Calculation | Status |
|-----------|---------|---------------------|-------------------|--------|
| E_MJ | Miyazawa-Jernigan interaction | ✓ Correct | ✓ Correct | ✓ Match |
| E1 | One position per residue | ✓ Correct | ✓ Correct | ✓ Match |
| E2 | Self-avoidance constraint | ✓ Correct | ✓ Correct | ✓ Match |
| E3 | Chain connectivity | ✓ Correct | ✓ Correct | ✓ Match |

---

## Detailed Mathematical Verification

### 1. E_MJ: Miyazawa-Jernigan Interaction Energy

**Formula:** `E_MJ = -Σ_{|i-j|>1} C(a_i, a_j) Σ_{<n,m>} b_{i,n} b_{j,m}`

**Implementation Analysis:**
- ✓ Correctly iterates over non-consecutive residue pairs (j >= i+2)
- ✓ Properly applies negative sign from formula (line 45 in qubo_generation.py)
- ✓ Only considers adjacent positions via adjacency matrix
- ✓ Upper triangular matrix storage (lines 55-58)

**Validation:** Tested with H-H favorable interactions (C=+1) and C-C repulsive interactions (C=-1)

---

### 2. E1: One Site Per Amino Acid Constraint

**Formula:** `E1 = Σ_i (Σ_n b_{i,n} - 1)²`

**Expanded Form:** `E1 = Σ_i [-Σ_n b_{i,n} + 2Σ_{n<m} b_{i,n}b_{i,m} + 1]`

**Implementation Analysis:**
- ✓ Diagonal terms: -1 for each b_{i,n} (lines 143-146)
- ✓ Off-diagonal terms: +2 for pairs within same residue (lines 149-154)
- ✓ Constant: +1 per residue (line 157)
- ✓ Correct expansion of (Σb - 1)²

**Key Insight:** When a residue occupies exactly one site, one b_{i,n}=1 and all others are 0, so:
- Linear terms: -1(1) = -1
- Quadratic terms: 2(0)(0) = 0 
- Constant: +1
- **Total: E1 = 0** ✓ (satisfied constraint)

---

### 3. E2: Self-Avoidance Constraint

**Formula:** `E2 = (1/2) Σ_n Σ_{i≠j} b_{i,n} b_{j,n}`

**Simplified:** `E2 = Σ_n Σ_{i<j} b_{i,n} b_{j,n}`

**Implementation Analysis:**
- ✓ Iterates over all positions n (line 237)
- ✓ For each position, considers all residue pairs i<j (lines 239-240)
- ✓ Coefficient is 1 in Q matrix (line 245)

---

### 4. E3: Chain Connectivity Constraint

**Formula:** `E3 = Σ_{i=0}^{N-2} Σ_n b_{i,n} Σ_{m: non-adj(n,m)} b_{i+1,m}`

**Implementation Analysis:**
- ✓ Constructs non-adjacency matrix: `non_adj = 1 - adj - I` (line 319)
- ✓ Correctly removes diagonal (positions can't be non-adjacent to themselves)
- ✓ Iterates over consecutive residue pairs (line 322)
- ✓ Penalizes non-adjacent position assignments (lines 324-337)

**Validation:** When consecutive residues are at adjacent positions, no penalty. When they're not adjacent, E3 > 0.

---

## Code Quality Assessment

### Strengths

1. **Modularity:** Excellent separation of concerns
   - QUBO construction (qubo_generation.py)
   - Energy evaluation (calc_mods.py)
   - Integration (protein_folding_qubo.py)

2. **Documentation:** 
   - Clear docstrings
   - Mathematical formulas in comments
   - Comprehensive notebook with explanations

3. **Validation:** 
   - Dual calculation methods (QUBO vs direct)
   - Multiple test cases covering edge cases
   - Explicit constraint checking

4. **Correctness:**
   - All mathematical transformations are accurate
   - Proper handling of binary variable properties (b² = b)
   - Correct QUBO matrix structure (upper triangular)

### Minor Issues and Recommendations

1. **E2 Calculation in calc_mods.py (lines 41-49):**
   - **Current:** Manually counts pairs and multiplies by 0.5
   - **Why it works:** This correctly implements (1/2) Σ_{i≠j}
   - **Note:** The formula counts ordered pairs then divides by 2
   - **Status:** Correct, but comment could clarify the double-counting logic

2. **Consistency:**
   - QUBO uses coefficient 1 with i<j ordering
   - Direct calculation uses 0.5 with explicit pair counting
   - **Status:** Both approaches are mathematically equivalent ✓

---

## Mathematical Consistency Check

### QUBO Form vs Direct Calculation

For each component, I verified that:

```
b^T Q b + constant = Direct_Energy_Calculation(b)
```

**Test 1 Results (Valid Configuration):**
- E_MJ: -1.00 (QUBO) = -1.00 (Direct) ✓
- E1: 0.00 (QUBO) = 0.00 (Direct) ✓
- E2: 0.00 (QUBO) = 0.00 (Direct) ✓
- E3: 0.00 (QUBO) = 0.00 (Direct) ✓
- Total: -1.00 (QUBO) = -1.00 (Direct) ✓

**Test 2 Results (E3 Violation):**
- E_MJ: -2.00 (QUBO) = -2.00 (Direct) ✓
- E1: 0.00 (QUBO) = 0.00 (Direct) ✓
- E2: 0.00 (QUBO) = 0.00 (Direct) ✓
- E3: 3.00 (QUBO) = 3.00 (Direct) ✓
- Total: 1.00 (QUBO) = 1.00 (Direct) ✓

**Test 3 Results (Multiple Violations):**
- Total: 8.00 (QUBO) = 8.00 (Direct) ✓
- E1 = 3 (violates one-position-per-residue)
- E2 = 4 (violates self-avoidance)
- E3 = 2 (violates connectivity)

---

## Constraint Behavior Verification

### E1 Constraint (One Position Per Residue)

- ✓ E1 = 0 when each residue at exactly one position
- ✓ E1 > 0 when residue at multiple positions
- ✓ E1 > 0 when residue at no position

### E2 Constraint (Self-Avoidance)

- ✓ E2 = 0 when no two residues share a position
- ✓ E2 > 0 when multiple residues at same position
- ✓ Correctly counts overlapping pairs

### E3 Constraint (Chain Connectivity)

- ✓ E3 = 0 when consecutive residues at adjacent positions
- ✓ E3 > 0 when consecutive residues not adjacent
- ✓ Non-adjacency matrix properly constructed

---

## Scalability Analysis

From protein_folding_qubo.py (lines 783-820):

| Residues | Positions | Variables | Matrix Size | Non-zero Terms | Sparsity |
|----------|-----------|-----------|-------------|----------------|----------|
| 4 | 4 | 16 | 256 | 84 | 67.2% |
| 6 | 6 | 36 | 1,296 | 336 | 74.1% |
| 8 | 8 | 64 | 4,096 | 768 | 81.2% |
| 10 | 10 | 100 | 10,000 | 1,380 | 86.2% |

**Observations:**
- ✓ Complexity scales as O(N²M²) as expected
- ✓ Matrices remain sparse (>65% zeros)
- ✓ Suitable for quantum annealing hardware

---

## Edge Cases Tested

1. ✓ Valid configuration with minimum energy
2. ✓ Single constraint violation (E3 only)
3. ✓ Multiple constraint violations (E1, E2, E3)
4. ✓ Different chain lengths (4, 6 residues)
5. ✓ Different amino acid compositions (H, P, C)
6. ✓ Various lattice topologies

---

## Potential Issues and Clarifications

### Issue 1: E2 Coefficient Discrepancy (Resolved)

**Location:** qubo_generation.py line 246, calc_mods.py line 49

**Analysis:**
- QUBO matrix uses coefficient 1 with i<j ordering
- Direct calculation uses 0.5 with full i≠j summation
- Both are mathematically equivalent

**Mathematics:**
```
E2 = (1/2) Σ_n Σ_{i≠j} b_{i,n} b_{j,n}
   = (1/2) Σ_n Σ_{i<j} 2·b_{i,n} b_{j,n}    [by symmetry]
   = Σ_n Σ_{i<j} b_{i,n} b_{j,n}            [QUBO form]
```

**Conclusion:** Implementation is correct; comment could be clearer.

---

## Recommendations

### Critical (None)
No critical issues found.

### Minor Improvements

3. **Add validation check** for adjacency matrix symmetry
   ```python
   def build_E3(chain, adj):
       # Verify adjacency matrix is symmetric
       assert np.allclose(adj, adj.T), "Adjacency matrix must be symmetric"
       ...
   ```

---

## Conclusion

### Overall Assessment: ✓ EXCELLENT

Your QUBO formulation for protein folding is **mathematically rigorous, correctly implemented, and thoroughly tested**. The code demonstrates:

1. ✓ Correct mathematical transformations from continuous to QUBO form
2. ✓ Proper constraint encoding as penalty terms
3. ✓ Accurate energy calculations via both methods
4. ✓ Comprehensive validation across multiple test cases
5. ✓ Clear documentation and code structure
6. ✓ Scalability analysis and complexity characterization

### Key Strengths

- **Mathematical Correctness:** All formulas correctly implemented
- **Validation:** Dual calculation methods provide strong verification
- **Code Quality:** Well-structured, documented, and modular
- **Completeness:** Handles all constraint types and interactions

---

## Test Evidence

All verification tests passed with exact numerical agreement between QUBO and direct calculation methods:

```
Test 1 (Valid config):     ✓ PASS
Test 2 (Invalid config):   ✓ PASS  
Test 3 (Multiple violat.): ✓ PASS

✓ ALL TESTS PASSED - IMPLEMENTATION IS CORRECT
```

---

**Final Verdict:** Your implementation is production-ready and mathematically sound. It correctly formulates the protein folding problem as a QUBO suitable for quantum annealing or classical optimization methods.