"""
Test suite for QUBO construction functions.
Validates that the QUBO matrices and polynomials are built correctly
by comparing computed energies against expected values.
"""

import numpy as np # type: ignore
from qubo_generation import (
    build_E_MJ,
    build_E1,
    build_E2,
    build_E3,
    print_E_MJ_details,
    print_E1_details,
    print_E2_details,
    print_E3_details,
)


# ============================================================================
# Testing and Validation
# ============================================================================

def evaluate_qubo_energy(bitstring, Q, constant):
    """
    Evaluate QUBO energy for a given bitstring.
    E = constant + sum_i Q[i,i]*x_i + sum_{i<j} 2*Q[i,j]*x_i*x_j
    """
    x = np.array([int(b) for b in bitstring.replace(' ', '')])
    
    energy = constant
    
    # Diagonal terms
    for i in range(len(x)):
        energy += Q[i, i] * x[i]
    
    # Off-diagonal terms (factor of 2 for QUBO convention)
    for i in range(len(x)):
        for j in range(i + 1, len(x)):
            energy += Q[i, j] * x[i] * x[j]
    
    return energy


def test_qubo_components(chain, adj, interaction_matrix, bitstring, expected_results, test_name):
    """
    Test QUBO construction by evaluating individual energy components.
    
    Args:
        chain: list of residue types
        adj: adjacency matrix
        interaction_matrix: interaction energies
        bitstring: binary string (spaces allowed)
        expected_results: tuple (E_MJ, E1, E2, E3)
        test_name: string identifier for this test
    """
    print("="*80)
    print(f"TEST: {test_name}")
    print("="*80)
    print(f"Chain: {chain}")
    print(f"Bitstring: {bitstring}")
    print(f"Expected: E_MJ={expected_results[0]}, E1={expected_results[1]}, "
          f"E2={expected_results[2]}, E3={expected_results[3]}")
    print()
    
    num_positions = adj.shape[0]
    
    # Build all components
    Q_MJ, _, const_MJ = build_E_MJ(chain, adj, interaction_matrix)
    Q_E1, _, const_E1 = build_E1(chain, num_positions)
    Q_E2, _, const_E2 = build_E2(chain, num_positions)
    Q_E3, _, const_E3 = build_E3(chain, adj)
    
    # Evaluate each component
    E_MJ_computed = evaluate_qubo_energy(bitstring, Q_MJ, const_MJ)
    E1_computed = evaluate_qubo_energy(bitstring, Q_E1, const_E1)
    E2_computed = evaluate_qubo_energy(bitstring, Q_E2, const_E2)
    E3_computed = evaluate_qubo_energy(bitstring, Q_E3, const_E3)
    
    print("RESULTS:")
    print(f"  E_MJ: computed = {E_MJ_computed:6.1f}, expected = {expected_results[0]:6.1f}  {'✓' if abs(E_MJ_computed - expected_results[0]) < 0.01 else '✗'}")
    print(f"  E1:   computed = {E1_computed:6.1f}, expected = {expected_results[1]:6.1f}  {'✓' if abs(E1_computed - expected_results[1]) < 0.01 else '✗'}")
    print(f"  E2:   computed = {E2_computed:6.1f}, expected = {expected_results[2]:6.1f}  {'✓' if abs(E2_computed - expected_results[2]) < 0.01 else '✗'}")
    print(f"  E3:   computed = {E3_computed:6.1f}, expected = {expected_results[3]:6.1f}  {'✓' if abs(E3_computed - expected_results[3]) < 0.01 else '✗'}")
    print()
    
    total_computed = E_MJ_computed + E1_computed + E2_computed + E3_computed
    total_expected = sum(expected_results)
    
    print(f"TOTAL ENERGY:")
    print(f"  Computed: {total_computed:.1f}")
    print(f"  Expected: {total_expected:.1f}")
    print(f"  Status: {'✓ PASS' if abs(total_computed - total_expected) < 0.01 else '✗ FAIL'}")
    print()
    
    return {
        'E_MJ': E_MJ_computed,
        'E1': E1_computed,
        'E2': E2_computed,
        'E3': E3_computed,
        'total': total_computed
    }


def run_all_tests():
    """Run all three test cases."""
    
    # Common interaction matrix
    C = {
        ('H', 'H'): 1, ('C', 'C'): -1,
        ('H', 'C'): 0, ('C', 'H'): 0,
        ('H', 'P'): 0, ('P', 'H'): 0,
        ('C', 'P'): 0, ('P', 'C'): 0,
        ('P', 'P'): 0
    }
    
    print("\n" + "="*80)
    print("RUNNING ALL QUBO TESTS")
    print("="*80)
    print()
    
    # ========================================================================
    # Test 1: Valid configuration
    # ========================================================================
    chain1 = ['H', 'P', 'C', 'H']
    adj1 = np.array([[0, 1, 1, 0],
                     [1, 0, 0, 1],
                     [1, 0, 0, 1],
                     [0, 1, 1, 0]])
    bitstring1 = '0010 1000 0100 0001'
    expected1 = (-1, 0, 0, 0)  # E_MJ, E1, E2, E3
    
    results1 = test_qubo_components(chain1, adj1, C, bitstring1, expected1, 
                                     "Test 1: Valid Configuration (HPCH)")
    
    # ========================================================================
    # Test 2: Multiple violations
    # ========================================================================
    chain2 = ['H', 'C', 'H', 'P']
    adj2 = np.array([[0, 1, 1, 0],
                     [1, 0, 0, 1],
                     [1, 0, 0, 1],
                     [0, 1, 1, 0]])
    bitstring2 = '1100 0101 0100 1010'
    expected2 = (-1, 3, 4, 2)  # E_MJ, E1, E2, E3
    
    results2 = test_qubo_components(chain2, adj2, C, bitstring2, expected2,
                                     "Test 2: Multiple Violations (HPCH)")
    
    # ========================================================================
    # Test 3: Larger chain (HHCHPC)
    # ========================================================================
    chain3 = ['H', 'H', 'C', 'H', 'P', 'C']
    adj3 = np.array([[0, 1, 1, 0, 0, 0],
                     [1, 0, 0, 1, 0, 0],
                     [1, 0, 0, 1, 1, 0],
                     [0, 1, 1, 0, 0, 1],
                     [0, 0, 1, 0, 0, 1],
                     [0, 0, 0, 1, 1, 0]])
    bitstring3 = '100000 000010 010000 001000 000100 000001'
    expected3 = (-2, 0, 0, 3)  # E_MJ, E1, E2, E3
    
    # First show the construction for this larger example
    print("="*80)
    print("DETAILED CONSTRUCTION: Test 3 (HHCHPC)")
    print("="*80)
    print(f"Chain: {chain3}")
    print(f"Number of residues: {len(chain3)}")
    print(f"Number of positions: {adj3.shape[0]}")
    print(f"Total bits: {len(chain3) * adj3.shape[0]}")
    print()
    print("Position Adjacency Matrix:")
    print(adj3)
    print()
    
    Q_MJ3, poly_MJ3, const_MJ3 = print_E_MJ_details(chain3, adj3, C)
    Q_E1_3, poly_E1_3, const_E1_3 = print_E1_details(chain3, adj3.shape[0])
    Q_E2_3, poly_E2_3, const_E2_3 = print_E2_details(chain3, adj3.shape[0])
    Q_E3_3, poly_E3_3, const_E3_3 = print_E3_details(chain3, adj3)
    
    print("="*80)
    print("SUMMARY FOR TEST 3")
    print("="*80)
    print(f"E_MJ: {len(poly_MJ3)} terms, constant = {const_MJ3}")
    print(f"E1:   {len(poly_E1_3)} terms, constant = {const_E1_3}")
    print(f"E2:   {len(poly_E2_3)} terms, constant = {const_E2_3}")
    print(f"E3:   {len(poly_E3_3)} terms, constant = {const_E3_3}")
    print()
    
    # Now test with the bitstring
    results3 = test_qubo_components(chain3, adj3, C, bitstring3, expected3,
                                     "Test 3: Larger Chain (HHCHPC)")
    
    # ========================================================================
    # Summary
    # ========================================================================
    print("="*80)
    print("FINAL SUMMARY")
    print("="*80)
    all_passed = (
        abs(results1['total'] - sum(expected1)) < 0.01 and
        abs(results2['total'] - sum(expected2)) < 0.01 and
        abs(results3['total'] - sum(expected3)) < 0.01
    )
    
    if all_passed:
        print("✓ ALL TESTS PASSED!")
    else:
        print("✗ SOME TESTS FAILED")
    print()


if __name__ == "__main__":
    run_all_tests()