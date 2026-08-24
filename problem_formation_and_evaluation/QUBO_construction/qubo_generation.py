import numpy as np # type: ignore
from itertools import combinations

def bit_index(residue_idx, position_idx, num_positions):
    """Convert (residue, position) to bit index."""
    return residue_idx * num_positions + position_idx

def decode_bit_index(bit_idx, num_positions):
    """Convert bit index to (residue, position)."""
    residue_idx = bit_idx // num_positions
    position_idx = bit_idx % num_positions
    return residue_idx, position_idx

# ============================================================================
# E_MJ: Miyazawa-Jernigan Interaction Energy
# ============================================================================

def build_E_MJ(chain, adj, interaction_matrix):
    """
    Build Q_MJ matrix and polynomial for Miyazawa-Jernigan interactions.
    
    E_MJ = -sum_{|i-j|>1} C(a_i, a_j) * sum_{<n,m>} b_{i,n} * b_{j,m}
    
    Returns:
        Q_MJ: numpy array (num_bits x num_bits)
        polynomial: list of (coefficient, bit_i, bit_j) tuples
        constant: constant offset
    """
    num_residues = len(chain)
    num_positions = adj.shape[0]
    num_bits = num_residues * num_positions
    
    Q_MJ = np.zeros((num_bits, num_bits))
    polynomial = []
    
    # Iterate over non-consecutive residue pairs
    for i in range(num_residues):
        for j in range(i + 2, num_residues):  # |i-j| > 1
            # Get interaction energy
            interaction = interaction_matrix.get((chain[i], chain[j]), 0)
            if interaction == 0:
                continue
            
            # Apply negative sign from E_MJ formula
            coeff = -interaction
            
            # For each pair of adjacent positions
            for n in range(num_positions):
                for m in range(num_positions):
                    if adj[n, m] == 1:  # Positions n and m are adjacent
                        bit_i = bit_index(i, n, num_positions)
                        bit_j = bit_index(j, m, num_positions)
                        
                        # Add to Q matrix (upper triangular)
                        if bit_i < bit_j:
                            Q_MJ[bit_i, bit_j] += coeff
                        else:
                            Q_MJ[bit_j, bit_i] += coeff
                        
                        # Add to polynomial
                        polynomial.append((coeff, bit_i, bit_j))
    
    constant = 0
    return Q_MJ, polynomial, constant


def print_E_MJ_details(chain, adj, interaction_matrix):
    """Print detailed breakdown of E_MJ construction."""
    num_residues = len(chain)
    num_positions = adj.shape[0]
    
    print("="*70)
    print("E_MJ: Miyazawa-Jernigan Interaction Energy")
    print("="*70)
    print(f"Formula: E_MJ = -sum_{{|i-j|>1}} C(a_i, a_j) * sum_{{<n,m>}} b_{{i,n}} * b_{{j,m}}")
    print()
    
    # Find non-consecutive pairs with non-zero interactions
    print("Non-consecutive residue pairs:")
    for i in range(num_residues):
        for j in range(i + 2, num_residues):
            interaction = interaction_matrix.get((chain[i], chain[j]), 0)
            print(f"  Pair ({i}, {j}): {chain[i]} and {chain[j]} -> C = {interaction}")
    print()
    
    # Build matrix
    Q_MJ, polynomial, constant = build_E_MJ(chain, adj, interaction_matrix)
    
    print(f"Q_MJ Matrix ({Q_MJ.shape[0]}x{Q_MJ.shape[1]}):")
    print(Q_MJ)
    print()
    
    print(f"Non-zero entries in Q_MJ:")
    for i in range(Q_MJ.shape[0]):
        for j in range(i, Q_MJ.shape[1]):
            if Q_MJ[i, j] != 0:
                res_i, pos_i = decode_bit_index(i, num_positions)
                res_j, pos_j = decode_bit_index(j, num_positions)
                print(f"  Q[{i:2d}, {j:2d}] = {Q_MJ[i,j]:+.1f}  (b_{{{res_i},{pos_i}}} * b_{{{res_j},{pos_j}}})")
    print()
    
    print(f"Polynomial (showing first 10 terms):")
    print("E_MJ = ", end="")
    for idx, (coeff, bit_i, bit_j) in enumerate(polynomial[:10]):
        res_i, pos_i = decode_bit_index(bit_i, num_positions)
        res_j, pos_j = decode_bit_index(bit_j, num_positions)
        if idx > 0:
            print(f"       {coeff:+.1f} * b_{{{res_i},{pos_i}}} * b_{{{res_j},{pos_j}}}")
        else:
            print(f"{coeff:+.1f} * b_{{{res_i},{pos_i}}} * b_{{{res_j},{pos_j}}}")
    if len(polynomial) > 10:
        print(f"       ... ({len(polynomial) - 10} more terms)")
    print()
    
    return Q_MJ, polynomial, constant


# ============================================================================
# E1: One Site Per Amino Acid
# ============================================================================

def build_E1(chain, num_positions):
    """
    Build Q_E1 matrix and polynomial.
    
    E1 = sum_i (sum_n b_{i,n} - 1)^2
       = sum_i [-sum_n b_{i,n} + 2*sum_{n<m} b_{i,n}*b_{i,m} + 1]
    
    Returns:
        Q_E1: numpy array (num_bits x num_bits)
        polynomial: list of (coefficient, bit_i, bit_j) tuples
        constant: constant offset
    """
    num_residues = len(chain)
    num_bits = num_residues * num_positions
    
    Q_E1 = np.zeros((num_bits, num_bits))
    polynomial = []
    constant = 0
    
    for i in range(num_residues):
        # Diagonal terms: -1 for each position
        for n in range(num_positions):
            bit_idx = bit_index(i, n, num_positions)
            Q_E1[bit_idx, bit_idx] += -1
            polynomial.append((-1, bit_idx, bit_idx))
        
        # Off-diagonal terms: +2 for each pair of positions within same residue
        for n in range(num_positions):
            for m in range(n + 1, num_positions):
                bit_n = bit_index(i, n, num_positions)
                bit_m = bit_index(i, m, num_positions)
                Q_E1[bit_n, bit_m] += 2
                polynomial.append((2, bit_n, bit_m))
        
        # Constant: +1 per residue
        constant += 1
    
    return Q_E1, polynomial, constant


def print_E1_details(chain, num_positions):
    """Print detailed breakdown of E1 construction."""
    num_residues = len(chain)
    
    print("="*70)
    print("E1: One Site Per Amino Acid")
    print("="*70)
    print(f"Formula: E1 = sum_i (sum_n b_{{i,n}} - 1)^2")
    print()
    
    Q_E1, polynomial, constant = build_E1(chain, num_positions)
    
    print(f"Q_E1 Matrix ({Q_E1.shape[0]}x{Q_E1.shape[1]}):")
    print(Q_E1)
    print()
    
    print(f"Diagonal entries (linear terms):")
    for i in range(Q_E1.shape[0]):
        if Q_E1[i, i] != 0:
            res_i, pos_i = decode_bit_index(i, num_positions)
            print(f"  Q[{i:2d}, {i:2d}] = {Q_E1[i,i]:+.1f}  (b_{{{res_i},{pos_i}}})")
    print()
    
    print(f"Off-diagonal entries (quadratic terms, showing first 10):")
    count = 0
    for i in range(Q_E1.shape[0]):
        for j in range(i+1, Q_E1.shape[1]):
            if Q_E1[i, j] != 0:
                res_i, pos_i = decode_bit_index(i, num_positions)
                res_j, pos_j = decode_bit_index(j, num_positions)
                print(f"  Q[{i:2d}, {j:2d}] = {Q_E1[i,j]:+.1f}  (b_{{{res_i},{pos_i}}} * b_{{{res_j},{pos_j}}})")
                count += 1
                if count >= 10:
                    break
        if count >= 10:
            break
    
    # Count total off-diagonal
    total_off_diag = np.sum(Q_E1 != 0) - num_residues * num_positions
    if total_off_diag > 10:
        print(f"  ... ({total_off_diag - 10} more off-diagonal terms)")
    print()
    
    print(f"Constant: {constant}")
    print()
    
    return Q_E1, polynomial, constant


# ============================================================================
# E2: Self-Avoidance (One Residue Per Position)
# ============================================================================

def build_E2(chain, num_positions):
    """
    Build Q_E2 matrix and polynomial.
    
    E2 = (1/2) * sum_n sum_{i!=j} b_{i,n} * b_{j,n}
       = sum_n sum_{i<j} b_{i,n} * b_{j,n}

    Implementation: We count site_occupancy*(site_occupancy-1) which gives
    the number of ordered pairs, then multiply by 0.5 to get the formula result.
    
    Returns:
        Q_E2: numpy array (num_bits x num_bits)
        polynomial: list of (coefficient, bit_i, bit_j) tuples
        constant: constant offset

        note: since b^2 = b for binary variables, any permutations wont show up at separate or multiplied terms
        Ex. since we have b_{0,1} * b_{1,1}, we'd technically also have b_{1,1} * b_{0,1}, but we'd just combine them, which evaluates to the same thing.
    """
    num_residues = len(chain)
    num_bits = num_residues * num_positions
    
    Q_E2 = np.zeros((num_bits, num_bits))
    polynomial = []
    
    # For each position
    for n in range(num_positions):
        # For each pair of residues at that position
        for i in range(num_residues):
            for j in range(i + 1, num_residues):
                bit_i = bit_index(i, n, num_positions)
                bit_j = bit_index(j, n, num_positions)
                
                # Add 1 to Q matrix (using i<j ordering eliminates need for 1/2 factor)
                Q_E2[bit_i, bit_j] += 1

                # Polynomial stores 0.5 for documentation (standard formula form)
                polynomial.append((0.5, bit_i, bit_j))
    
    constant = 0
    return Q_E2, polynomial, constant


def print_E2_details(chain, num_positions):
    """Print detailed breakdown of E2 construction."""
    num_residues = len(chain)
    
    print("="*70)
    print("E2: Self-Avoidance (One Residue Per Position)")
    print("="*70)
    print(f"Formula: E2 = (1/2) * sum_n sum_{{i!=j}} b_{{i,n}} * b_{{j,n}}")
    print()
    
    Q_E2, polynomial, constant = build_E2(chain, num_positions)
    
    print(f"Q_E2 Matrix ({Q_E2.shape[0]}x{Q_E2.shape[1]}):")
    print(Q_E2)
    print()
    
    print(f"Off-diagonal entries (showing first 10):")
    count = 0
    for i in range(Q_E2.shape[0]):
        for j in range(i+1, Q_E2.shape[1]):
            if Q_E2[i, j] != 0:
                res_i, pos_i = decode_bit_index(i, num_positions)
                res_j, pos_j = decode_bit_index(j, num_positions)
                print(f"  Q[{i:2d}, {j:2d}] = {Q_E2[i,j]:+.1f}  (b_{{{res_i},{pos_i}}} * b_{{{res_j},{pos_j}}})")
                count += 1
                if count >= 10:
                    break
        if count >= 10:
            break
    
    total_off_diag = np.sum(Q_E2 != 0)
    if total_off_diag > 10:
        print(f"  ... ({total_off_diag - 10} more off-diagonal terms)")
    print()
    
    print(f"Note: All diagonal entries are 0 (E2 is purely quadratic)")
    print(f"Constant: {constant}")
    print()
    
    return Q_E2, polynomial, constant


# ============================================================================
# E3: Chain Connectivity
# ============================================================================

def build_E3(chain, adj):
    """
    Build Q_E3 matrix and polynomial.
    
    E3 = sum_{i=1}^{N-1} sum_n b_{i,n} sum_{||m-n||>1} b_{i+1,m}
    
    This penalizes consecutive residues that are NOT at adjacent positions.
    
    Returns:
        Q_E3: numpy array (num_bits x num_bits)
        polynomial: list of (coefficient, bit_i, bit_j) tuples
        constant: constant offset
    """
    num_residues = len(chain)
    num_positions = adj.shape[0]
    num_bits = num_residues * num_positions
    
    Q_E3 = np.zeros((num_bits, num_bits))
    polynomial = []
    
    # Build non-adjacency matrix (excluding diagonal)
    non_adj = 1 - adj - np.eye(num_positions)
    
    # For each consecutive residue pair
    for i in range(num_residues - 1):
        # For each pair of non-adjacent positions (n, m)
        for n in range(num_positions):
            for m in range(num_positions):
                if non_adj[n, m] == 1:  # Positions n and m are NOT adjacent
                    bit_i = bit_index(i, n, num_positions)
                    bit_j = bit_index(i + 1, m, num_positions)
                    
                    # Add to Q matrix (upper triangular)
                    if bit_i < bit_j:
                        Q_E3[bit_i, bit_j] += 1
                    else:
                        Q_E3[bit_j, bit_i] += 1
                    
                    # Add to polynomial
                    polynomial.append((1, bit_i, bit_j))
    
    constant = 0
    return Q_E3, polynomial, constant


def print_E3_details(chain, adj):
    """Print detailed breakdown of E3 construction."""
    num_residues = len(chain)
    num_positions = adj.shape[0]
    
    print("="*70)
    print("E3: Chain Connectivity")
    print("="*70)
    print(f"Formula: E3 = sum_{{i=1}}^{{N-1}} sum_n b_{{i,n}} sum_{{||m-n||>1}} b_{{i+1,m}}")
    print()
    print("This penalizes consecutive residues at NON-adjacent positions.")
    print()
    
    # Show non-adjacent position pairs
    non_adj = 1 - adj - np.eye(num_positions)
    print("Non-adjacent position pairs:")
    for n in range(num_positions):
        for m in range(num_positions):
            if non_adj[n, m] == 1:
                print(f"  Positions {n} and {m} are NOT adjacent")
    print()
    
    print(f"Consecutive residue pairs: {num_residues - 1}")
    for i in range(num_residues - 1):
        print(f"  Pair ({i}, {i+1}): {chain[i]} -> {chain[i+1]}")
    print()
    
    Q_E3, polynomial, constant = build_E3(chain, adj)
    
    print(f"Q_E3 Matrix ({Q_E3.shape[0]}x{Q_E3.shape[1]}):")
    print(Q_E3)
    print()
    
    print(f"Non-zero entries in Q_E3:")
    for i in range(Q_E3.shape[0]):
        for j in range(i, Q_E3.shape[1]):
            if Q_E3[i, j] != 0:
                res_i, pos_i = decode_bit_index(i, num_positions)
                res_j, pos_j = decode_bit_index(j, num_positions)
                print(f"  Q[{i:2d}, {j:2d}] = {Q_E3[i,j]:+.1f}  (b_{{{res_i},{pos_i}}} * b_{{{res_j},{pos_j}}})  [res {res_i}->{res_j}, pos {pos_i}->{pos_j} non-adj]")
    print()
    
    print(f"Total terms: {len(polynomial)}")
    print(f"Constant: {constant}")
    print()
    
    return Q_E3, polynomial, constant


# ============================================================================
# Main demonstration
# ============================================================================

if __name__ == "__main__":
    # Example setup
    chain = ['H', 'P', 'C', 'H']
    adj = np.array([[0, 1, 1, 0],
                    [1, 0, 0, 1],
                    [1, 0, 0, 1],
                    [0, 1, 1, 0]])
    
    C = {
        ('H', 'H'): 1, ('C', 'C'): -1,
        ('H', 'C'): 0, ('C', 'H'): 0,
        ('H', 'P'): 0, ('P', 'H'): 0,
        ('C', 'P'): 0, ('P', 'C'): 0,
        ('P', 'P'): 0
    }
    
    num_positions = adj.shape[0]
    
    print("\n" + "="*70)
    print("QUBO FORMULATION FOR PROTEIN FOLDING")
    print("="*70)
    print(f"Chain: {chain}")
    print(f"Number of positions: {num_positions}")
    print(f"Total bits: {len(chain) * num_positions}")
    print()
    print("Position Adjacency Matrix:")
    print(adj)
    print()
    
    # Build and display E_MJ
    Q_MJ, poly_MJ, const_MJ = print_E_MJ_details(chain, adj, C)
    
    # Build and display E1
    Q_E1, poly_E1, const_E1 = print_E1_details(chain, num_positions)
    
    # Build and display E2
    Q_E2, poly_E2, const_E2 = print_E2_details(chain, num_positions)
    
    # Build and display E3
    Q_E3, poly_E3, const_E3 = print_E3_details(chain, adj)
    
    print("="*70)
    print("SUMMARY")
    print("="*70)
    print(f"E_MJ: {len(poly_MJ)} terms, constant = {const_MJ}")
    print(f"E1:   {len(poly_E1)} terms, constant = {const_E1}")
    print(f"E2:   {len(poly_E2)} terms, constant = {const_E2}")
    print(f"E3:   {len(poly_E3)} terms, constant = {const_E3}")
    print()
    
    # Combined totals
    total_terms = len(poly_MJ) + len(poly_E1) + len(poly_E2) + len(poly_E3)
    total_constant = const_MJ + const_E1 + const_E2 + const_E3
    print(f"TOTAL: {total_terms} terms, total constant = {total_constant}")