#!/usr/bin/env python3
"""
Comprehensive verification of QUBO protein folding implementation
"""

import numpy as np
import sys
import os

# Add the parent directory's subdirectories to path
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
qubo_construction_path = os.path.join(parent_dir, 'qubo_construction')
energy_calc_path = os.path.join(parent_dir, 'energy_calc')
sys.path.insert(0, qubo_construction_path)
sys.path.insert(0, energy_calc_path)

# Import the modules
from qubo_generation import ( #type: ignore
    build_E_MJ, build_E1, build_E2, build_E3,
    bit_index, decode_bit_index
)
from calc_mods import ( #type: ignore
    compute_E_MJ_debug, compute_E1, compute_E2, compute_E3,
    total_energy, is_valid_conformation, C
)

def evaluate_qubo_energy(bitstring, Q, offset):
    """Evaluate energy using QUBO matrix."""
    bitstring_clean = bitstring.replace(' ', '')
    b = np.array([int(bit) for bit in bitstring_clean])
    return b @ Q @ b + offset

def verify_component(name, chain, bitstring, adj, C, L_value, build_func, compute_func):
    """Verify a single energy component."""
    num_positions = adj.shape[0]
    
    # Build QUBO
    if name == "E_MJ":
        Q, poly, const = build_func(chain, adj, C)
    elif name == "E3":
        Q, poly, const = build_func(chain, adj)
    else:
        Q, poly, const = build_func(chain, num_positions)
    
    # Compute energy via QUBO
    E_qubo = evaluate_qubo_energy(bitstring, Q, const)
    
    # Compute energy directly
    N = len(chain)
    b = np.array([list(map(int, bitstring[i*N:(i+1)*N])) for i in range(N)])
    
    if name == "E_MJ":
        E_direct = compute_func(chain, b, adj, C, verbose=False)
    elif name == "E3":
        E_direct = compute_func(b, adj)
    else:
        E_direct = compute_func(b)
    
    match = abs(E_qubo - E_direct) < 0.01
    
    return {
        'name': name,
        'E_qubo': E_qubo,
        'E_direct': E_direct,
        'match': match,
        'num_terms': len(poly),
        'constant': const
    }

def run_comprehensive_tests():
    """Run comprehensive verification tests."""
    
    print("="*80)
    print("COMPREHENSIVE QUBO VERIFICATION")
    print("="*80)
    
    # Test 1: Simple valid configuration
    print("\n" + "="*80)
    print("TEST 1: Simple Valid Configuration (HPCH)")
    print("="*80)
    
    chain1 = ['H', 'P', 'C', 'H']
    bitstring1 = '0010100001000001'
    adj1 = np.array([[0, 1, 1, 0],
                     [1, 0, 0, 1],
                     [1, 0, 0, 1],
                     [0, 1, 1, 0]])
    
    results1 = []
    results1.append(verify_component("E_MJ", chain1, bitstring1, adj1, C, 1.0, build_E_MJ, compute_E_MJ_debug))
    results1.append(verify_component("E1", chain1, bitstring1, adj1, C, 1.0, build_E1, compute_E1))
    results1.append(verify_component("E2", chain1, bitstring1, adj1, C, 1.0, build_E2, compute_E2))
    results1.append(verify_component("E3", chain1, bitstring1, adj1, C, 1.0, build_E3, compute_E3))
    
    print(f"\nChain: {chain1}")
    print(f"Bitstring: {bitstring1}")
    print(f"\nComponent Verification:")
    print(f"{'Component':<10} {'QUBO Energy':<15} {'Direct Energy':<15} {'Match':<10} {'Terms':<10}")
    print("-"*70)
    
    all_match = True
    for r in results1:
        match_str = "✓" if r['match'] else "✗"
        print(f"{r['name']:<10} {r['E_qubo']:<15.2f} {r['E_direct']:<15.2f} {match_str:<10} {r['num_terms']:<10}")
        all_match = all_match and r['match']
    
    # Compute total energy
    total_E1, breakdown1 = total_energy(chain1, bitstring1, adj1, C, 1.0, 1.0, 1.0, verbose=False)
    Q_total1 = sum([build_E_MJ(chain1, adj1, C)[0],
                    build_E1(chain1, 4)[0],
                    build_E2(chain1, 4)[0],
                    build_E3(chain1, adj1)[0]])
    offset1 = sum([build_E_MJ(chain1, adj1, C)[2],
                   build_E1(chain1, 4)[2],
                   build_E2(chain1, 4)[2],
                   build_E3(chain1, adj1)[2]])
    E_qubo_total1 = evaluate_qubo_energy(bitstring1, Q_total1, offset1)
    
    print(f"\nTotal Energy:")
    print(f"  QUBO method:   {E_qubo_total1:.2f}")
    print(f"  Direct method: {total_E1:.2f}")
    print(f"  Match: {'✓' if abs(E_qubo_total1 - total_E1) < 0.01 else '✗'}")
    print(f"\nValid conformation: {is_valid_conformation(breakdown1)}")
    
    test1_pass = all_match and abs(E_qubo_total1 - total_E1) < 0.01
    
    # Test 2: Invalid configuration with constraint violations
    print("\n" + "="*80)
    print("TEST 2: Invalid Configuration with Violations (HHCHPC)")
    print("="*80)
    
    chain2 = ['H', 'H', 'C', 'H', 'P', 'C']
    bitstring2 = '100000000010010000001000000100000001'
    adj2 = np.array([[0, 1, 1, 0, 0, 0],
                     [1, 0, 0, 1, 0, 0],
                     [1, 0, 0, 1, 1, 0],
                     [0, 1, 1, 0, 0, 1],
                     [0, 0, 1, 0, 0, 1],
                     [0, 0, 0, 1, 1, 0]])
    
    results2 = []
    results2.append(verify_component("E_MJ", chain2, bitstring2, adj2, C, 1.0, build_E_MJ, compute_E_MJ_debug))
    results2.append(verify_component("E1", chain2, bitstring2, adj2, C, 1.0, build_E1, compute_E1))
    results2.append(verify_component("E2", chain2, bitstring2, adj2, C, 1.0, build_E2, compute_E2))
    results2.append(verify_component("E3", chain2, bitstring2, adj2, C, 1.0, build_E3, compute_E3))
    
    print(f"\nChain: {chain2}")
    print(f"Bitstring: {bitstring2}")
    print(f"\nComponent Verification:")
    print(f"{'Component':<10} {'QUBO Energy':<15} {'Direct Energy':<15} {'Match':<10} {'Terms':<10}")
    print("-"*70)
    
    all_match2 = True
    for r in results2:
        match_str = "✓" if r['match'] else "✗"
        print(f"{r['name']:<10} {r['E_qubo']:<15.2f} {r['E_direct']:<15.2f} {match_str:<10} {r['num_terms']:<10}")
        all_match2 = all_match2 and r['match']
    
    total_E2, breakdown2 = total_energy(chain2, bitstring2, adj2, C, 1.0, 1.0, 1.0, verbose=False)
    Q_total2 = sum([build_E_MJ(chain2, adj2, C)[0],
                    build_E1(chain2, 6)[0],
                    build_E2(chain2, 6)[0],
                    build_E3(chain2, adj2)[0]])
    offset2 = sum([build_E_MJ(chain2, adj2, C)[2],
                   build_E1(chain2, 6)[2],
                   build_E2(chain2, 6)[2],
                   build_E3(chain2, adj2)[2]])
    E_qubo_total2 = evaluate_qubo_energy(bitstring2, Q_total2, offset2)
    
    print(f"\nTotal Energy:")
    print(f"  QUBO method:   {E_qubo_total2:.2f}")
    print(f"  Direct method: {total_E2:.2f}")
    print(f"  Match: {'✓' if abs(E_qubo_total2 - total_E2) < 0.01 else '✗'}")
    print(f"\nConstraint violations:")
    print(f"  E1 (one pos/residue): {breakdown2[1]} {'✗' if breakdown2[1] > 0 else '✓'}")
    print(f"  E2 (self-avoidance):  {breakdown2[2]} {'✗' if breakdown2[2] > 0 else '✓'}")
    print(f"  E3 (connectivity):    {breakdown2[3]} {'✗' if breakdown2[3] > 0 else '✓'}")
    print(f"\nValid conformation: {is_valid_conformation(breakdown2)}")
    
    test2_pass = all_match2 and abs(E_qubo_total2 - total_E2) < 0.01
    
    # Test 3: Multiple violations
    print("\n" + "="*80)
    print("TEST 3: Multiple Constraint Violations (HCHP)")
    print("="*80)
    
    chain3 = ['H', 'C', 'H', 'P']
    bitstring3 = '1100010101001010'
    adj3 = np.array([[0, 1, 1, 0],
                     [1, 0, 0, 1],
                     [1, 0, 0, 1],
                     [0, 1, 1, 0]])
    
    total_E3, breakdown3 = total_energy(chain3, bitstring3, adj3, C, 1.0, 1.0, 1.0, verbose=False)
    Q_total3 = sum([build_E_MJ(chain3, adj3, C)[0],
                    build_E1(chain3, 4)[0],
                    build_E2(chain3, 4)[0],
                    build_E3(chain3, adj3)[0]])
    offset3 = sum([build_E_MJ(chain3, adj3, C)[2],
                   build_E1(chain3, 4)[2],
                   build_E2(chain3, 4)[2],
                   build_E3(chain3, adj3)[2]])
    E_qubo_total3 = evaluate_qubo_energy(bitstring3, Q_total3, offset3)
    
    print(f"\nChain: {chain3}")
    print(f"Bitstring: {bitstring3}")
    print(f"\nTotal Energy:")
    print(f"  QUBO method:   {E_qubo_total3:.2f}")
    print(f"  Direct method: {total_E3:.2f}")
    print(f"  Match: {'✓' if abs(E_qubo_total3 - total_E3) < 0.01 else '✗'}")
    print(f"\nConstraint violations:")
    print(f"  E1 (one pos/residue): {breakdown3[1]} {'✗' if breakdown3[1] > 0 else '✓'}")
    print(f"  E2 (self-avoidance):  {breakdown3[2]} {'✗' if breakdown3[2] > 0 else '✓'}")
    print(f"  E3 (connectivity):    {breakdown3[3]} {'✗' if breakdown3[3] > 0 else '✓'}")
    print(f"\nValid conformation: {is_valid_conformation(breakdown3)}")
    
    test3_pass = abs(E_qubo_total3 - total_E3) < 0.01
    
    # Final summary
    print("\n" + "="*80)
    print("VERIFICATION SUMMARY")
    print("="*80)
    print(f"\nTest 1 (Valid config):     {'✓ PASS' if test1_pass else '✗ FAIL'}")
    print(f"Test 2 (Invalid config):   {'✓ PASS' if test2_pass else '✗ FAIL'}")
    print(f"Test 3 (Multiple violat.): {'✓ PASS' if test3_pass else '✗ FAIL'}")
    
    if test1_pass and test2_pass and test3_pass:
        print("\n" + "="*80)
        print("✓ ALL TESTS PASSED - IMPLEMENTATION IS CORRECT")
        print("="*80)
        return True
    else:
        print("\n" + "="*80)
        print("✗ SOME TESTS FAILED - PLEASE REVIEW")
        print("="*80)
        return False

if __name__ == "__main__":
    success = run_comprehensive_tests()
    sys.exit(0 if success else 1)