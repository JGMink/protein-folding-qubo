"""
Regression test: the E_MJ objective must REWARD contacts, not penalize them.

Background
---------
Two sign conventions appear in the lattice-protein literature:

  (A)  E = -sum(C * contacts),  C = POSITIVE contact strength      (H-H = +1)
  (B)  E =  sum(J * contacts),  J = SIGNED Miyazawa-Jernigan energy (H-H < 0)

Both are correct in isolation. The bug is applying (A)'s negation to a table
that is already signed like (B): the double negation makes forming a contact
RAISE the energy, inverting the objective.

qubo_generation.py uses convention (A) consistently (C['H','H'] = +1, and
`coeff = -interaction`). These witnesses pin that down so a future edit that
reintroduces a double negation fails loudly.

Run directly (`python test_objective_sign.py`) or under pytest.
"""

import numpy as np

from qubo_generation import build_E_MJ, build_E1, build_E2, build_E3, bit_index

# Toy HP interaction matrix, convention (A): positive = favourable contact.
C = {
    ('H', 'H'): 1, ('C', 'C'): -1,
    ('H', 'C'): 0, ('C', 'H'): 0,
    ('H', 'P'): 0, ('P', 'H'): 0,
    ('C', 'P'): 0, ('P', 'C'): 0,
    ('P', 'P'): 0,
}


def _grid_2x3_adj():
    """2x3 lattice, sites laid out as:  0 1 2
                                        3 4 5
    Two sites are adjacent iff their Manhattan distance is 1.
    """
    coords = {0: (0, 0), 1: (1, 0), 2: (2, 0), 3: (0, 1), 4: (1, 1), 5: (2, 1)}
    adj = np.zeros((6, 6), dtype=int)
    for i in coords:
        for j in coords:
            if i == j:
                continue
            (x1, y1), (x2, y2) = coords[i], coords[j]
            if abs(x1 - x2) + abs(y1 - y2) == 1:
                adj[i, j] = 1
    return adj


def _energy_components(chain, sites, adj):
    """Place residue r at lattice site `sites[r]` and evaluate each QUBO block
    as x^T Q x + const (Q is upper-triangular, x binary)."""
    npos = adj.shape[0]
    nres = len(chain)
    x = np.zeros(nres * npos)
    for r, s in enumerate(sites):
        x[bit_index(r, s, npos)] = 1

    blocks = {
        'E_MJ': build_E_MJ(chain, adj, C),
        'E1': build_E1(chain, npos),
        'E2': build_E2(chain, npos),
        'E3': build_E3(chain, adj),
    }
    out = {name: float(x @ Q @ x + const) for name, (Q, _poly, const) in blocks.items()}
    out['total'] = sum(out.values())
    return out


def test_hpph_contact_fold_is_one_lower():
    """HPPH on a 2x3 grid.

    contact fold      sites [0, 1, 4, 3] -> residues 0 and 3 land on lattice-
                                            adjacent sites 0 and 3 (an H-H contact)
    non-contact fold  sites [0, 1, 2, 5] -> residues 0 and 3 not adjacent

    Both placements are valid self-avoiding connected chains, so E1 = E2 = E3 = 0
    and only E_MJ differs. The contact fold must be exactly 1.0 LOWER.
    """
    adj = _grid_2x3_adj()
    chain = list("HPPH")

    contact = _energy_components(chain, [0, 1, 4, 3], adj)
    non_contact = _energy_components(chain, [0, 1, 2, 5], adj)

    for comp in ('E1', 'E2', 'E3'):
        assert contact[comp] == 0.0, f"contact fold not valid: {comp}={contact[comp]}"
        assert non_contact[comp] == 0.0, f"non-contact fold not valid: {comp}={non_contact[comp]}"

    delta = contact['total'] - non_contact['total']
    assert delta == -1.0, (
        f"expected contact fold 1.0 lower, got delta={delta:+.1f}. "
        f"delta > 0 means the E_MJ objective is INVERTED (double negation)."
    )


def test_hph_ground_state_is_zero():
    """HPH: the lone non-consecutive pair (0, 2) can never be lattice-adjacent
    on a connected 3-residue chain (same parity => Manhattan distance >= 2), so
    the ground-state energy over valid placements is exactly 0. Any encoding
    that reports -1 for HPH has manufactured a geometrically impossible contact.
    """
    from itertools import permutations

    adj = _grid_2x3_adj()
    chain = list("HPH")

    best = None
    for sites in permutations(range(6), 3):
        ec = _energy_components(chain, list(sites), adj)
        if ec['E1'] == 0.0 and ec['E2'] == 0.0 and ec['E3'] == 0.0:
            if best is None or ec['total'] < best:
                best = ec['total']

    assert best is not None, "no valid HPH placement found on 2x3 grid"
    assert best == 0.0, f"HPH ground state should be 0.0, got {best:+.1f}"


if __name__ == "__main__":
    test_hpph_contact_fold_is_one_lower()
    test_hph_ground_state_is_zero()
    print("OK: E_MJ objective rewards contacts (HPPH delta = -1.0, HPH ground state = 0.0)")
