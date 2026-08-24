import numpy as np # type: ignore

# === Lagrange Penalties (Globals) ===
L1, L2, L3 = 1.0, 1.0, 1.0  # All set to 1.0 as requested

# === Your original interaction matrix ===
C = {
    ('H', 'H'): 1, ('C', 'C'): -1,
    ('H', 'C'): 0, ('C', 'H'): 0,
    ('H', 'P'): 0, ('P', 'H'): 0,
    ('C', 'P'): 0, ('P', 'C'): 0,
    ('P', 'P'): 0
}

def compute_E_MJ_debug(chain, b, adj, C, verbose=False):
    """
    Compute interaction energy with optional debug output
    """
    N = len(chain)
    E = 0.0
    contacts = []
    
    for i in range(N):
        for j in range(i+2, N):  # |i-j| > 1 means j >= i+2
            a_i, a_j = chain[i], chain[j]
            c_ij = C.get((a_i, a_j), 0.0)
            if c_ij != 0:
                for n in range(N):
                    for m in range(N):
                        if b[i, n] == 1 and b[j, m] == 1 and adj[n, m] == 1:
                            E += c_ij
                            contacts.append((i, j, a_i, a_j, n, m, c_ij))
                            if verbose:
                                print(f"  Contact: amino acids {i}({a_i}) and {j}({a_j}) at sites {n} and {m}, C={c_ij}")
    
    if verbose:
        print(f"  Total interaction sum: {E}")
        print(f"  E_MJ = -{E} = {-E}")
    
    return -E

def compute_E1(b):
    """
    Ensures each amino acid occupies exactly one site
    """
    result = np.sum((np.sum(b, axis=1) - 1) ** 2)
    return int(result)

def compute_E2(b):
    """
    Ensures each site has at most one amino acid
    E2 = (1/2) * sum_n sum_{i≠j} b_{i,n} b_{j,n}
    """
    N = b.shape[0]
    pair_count = 0
    for n in range(N):
        site_occupancy = np.sum(b[:, n])
        pair_count += site_occupancy * (site_occupancy - 1)
    # site_occupancy * (site_occupancy - 1) already counts ordered pairs
    # so we divide by 2 to get unordered pairs, then multiply by 1/2 from formula
    # Net effect: pair_count already has the factor of 2 built in
    return int(0.5 * pair_count)

def compute_E3(b, adj):
    """
    Ensures chain connectivity
    """
    N = b.shape[0]
    # Create non-adjacency matrix
    non_adj = 1 - adj - np.eye(N)
    
    E3 = 0.0
    for i in range(N - 1):
        E3 += b[i, :] @ non_adj @ b[i+1, :]
    
    return int(E3)

def total_energy(chain, bitstring, adj, C, L1=1.0, L2=1.0, L3=1.0, verbose=False):
    """
    Compute total QUBO energy with optional debug output
    """
    N = len(chain)
    b = np.array([list(map(int, bitstring[i*N:(i+1)*N])) for i in range(N)])
    
    if verbose:
        print("\nBit matrix b:")
        print(b)
        print("\nPositions:")
        for i in range(N):
            pos = np.where(b[i, :] == 1)[0]
            print(f"  Amino acid {i} ({chain[i]}): site {pos[0] if len(pos) > 0 else 'none'}")
        print("\nComputing E_MJ:")
    
    E_MJ = compute_E_MJ_debug(chain, b, adj, C, verbose)
    E_1  = compute_E1(b)
    E_2  = compute_E2(b)
    E_3  = compute_E3(b, adj)
    
    total = E_MJ + L1*E_1 + L2*E_2 + L3*E_3
    
    if verbose:
        print(f"\nEnergy components:")
        print(f"  E_MJ = {E_MJ}")
        print(f"  E_1  = {E_1}")
        print(f"  E_2  = {E_2}")
        print(f"  E_3  = {E_3}")
        print(f"\nTotal energy calculation:")
        print(f"  E = {E_MJ} + {L1}×{E_1} + {L2}×{E_2} + {L3}×{E_3}")
        print(f"  E = {E_MJ} + {L1*E_1} + {L2*E_2} + {L3*E_3}")
        print(f"  E = {total}")
    
    return total, (E_MJ, E_1, E_2, E_3)

def is_valid_conformation(breakdown):
    """Check if configuration satisfies all constraints"""
    _, E1, E2, E3 = breakdown
    return E1 == 0 and E2 == 0 and E3 == 0

if __name__ == "__main__":
    # Example usage 1 -> Simple and non-disruptive chain
    print("="*60)
    print("=== Example 1: Valid chain ===")
    print("="*60)
    chain = ['H', 'P', 'C', 'H']
    bitstring = '0010100001000001'
    adj = np.array([[0, 1, 1, 0],
                    [1, 0, 0, 1],
                    [1, 0, 0, 1],
                    [0, 1, 1, 0]])

    total_E, breakdown = total_energy(chain, bitstring, adj, C, L1, L2, L3, verbose=True)
    print(f"\nTotal Energy: {total_E}")
    print(f"Energy Breakdown: MJ={breakdown[0]}, E1={breakdown[1]}, E2={breakdown[2]}, E3={breakdown[3]}")
    print(f"Valid conformation: {is_valid_conformation(breakdown)}")

    # Example usage 2 -> just issues with chain connectivity, elsewise fine
    print("\n" + "="*60)
    print("=== Example 2: Chain connectivity violation ===")
    print("="*60)
    chain = ['H', 'H', 'C', 'H', 'P', 'C']
    bitstring = '100000000010010000001000000100000001'
    adj = np.array([[0, 1, 1, 0, 0, 0],
                    [1, 0, 0, 1, 0, 0],
                    [1, 0, 0, 1, 1, 0],
                    [0, 1, 1, 0, 0, 1],
                    [0, 0, 1, 0, 0, 1],
                    [0, 0, 0, 1, 1, 0]])

    total_E, breakdown = total_energy(chain, bitstring, adj, C, L1, L2, L3, verbose=True)
    print(f"\nTotal Energy: {total_E}")
    print(f"Energy Breakdown: MJ={breakdown[0]}, E1={breakdown[1]}, E2={breakdown[2]}, E3={breakdown[3]}")
    print(f"Valid conformation: {is_valid_conformation(breakdown)}")
    
    # Example usage 3 -> violates all 3 constraints, can't be represented physically
    print("\n" + "="*60)
    print("=== Example 3: Multiple violations ===")
    print("="*60)
    chain = ['H', 'C', 'H', 'P']
    bitstring = '1100010101001010'
    adj = np.array([[0, 1, 1, 0],
                    [1, 0, 0, 1],
                    [1, 0, 0, 1],
                    [0, 1, 1, 0]])

    total_E, breakdown = total_energy(chain, bitstring, adj, C, L1, L2, L3, verbose=True)
    print(f"\nTotal Energy: {total_E}")
    print(f"Energy Breakdown: MJ={breakdown[0]}, E1={breakdown[1]}, E2={breakdown[2]}, E3={breakdown[3]}")
    print(f"Valid conformation: {is_valid_conformation(breakdown)}")