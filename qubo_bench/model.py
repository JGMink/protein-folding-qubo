"""Sparse equivalent of the signed-contact notebook QUBO construction."""

from collections import defaultdict
import json
from pathlib import Path


DATA = json.loads((Path(__file__).parent / "contact_matrices.json").read_text(encoding="utf-8"))


def lattice(dims, order="row-major"):
    """Return coordinates and undirected edges in the notebook's site order."""
    if len(dims) not in (2, 3) or any(d < 1 for d in dims):
        raise ValueError("lattice dimensions must contain two or three positive integers")
    if len(dims) == 2:
        rows, cols = dims
        if order == "row-major":
            coords = [(r, c) for r in range(rows) for c in range(cols)]
        elif order == "column-major":
            coords = [(r, c) for c in range(cols) for r in range(rows)]
        else:
            raise ValueError("unsupported 2D site order")
    elif order == "x-fastest":
        x, y, z = dims
        coords = [(i, j, k) for k in range(z) for j in range(y) for i in range(x)]
    else:
        raise ValueError("3D lattices use x-fastest site order")
    edges = [(u, v) for u, a in enumerate(coords) for v in range(u + 1, len(coords))
             if sum(abs(p - q) for p, q in zip(a, coords[v])) == 1]
    return coords, edges


def build(sequence, dims, penalties=(3.0, 4.0, 4.0), model="signed20", order="row-major"):
    """Construct E_contact + L1*E1 + L2*E2 + L3*E3 as an upper QUBO.

    signed20 is the supplied signed notebook table. toy uses the tutorial's
    positive contact strengths after negation. No formulation is altered.
    """
    if not sequence:
        raise ValueError("sequence cannot be empty")
    if len(penalties) != 3 or any(p <= 0 for p in penalties):
        raise ValueError("three positive penalties are required")
    if model not in DATA:
        raise ValueError(f"unknown contact model: {model}")
    contact = DATA[model]
    alphabet = contact["alphabet"]
    indices = {aa: i for i, aa in enumerate(alphabet)}
    if any(aa not in indices for aa in sequence):
        raise ValueError("sequence contains residue absent from contact model")
    coords, edges = lattice(dims, order)
    n_res, n_sites = len(sequence), len(coords)
    l1, l2, l3 = map(float, penalties)
    q = defaultdict(float)

    def add(a, b, coefficient):
        q[min(a, b), max(a, b)] += coefficient

    for i in range(n_res):
        for j in range(i + 2, n_res):
            coefficient = contact["values"][indices[sequence[i]]][indices[sequence[j]]]
            if model == "toy":
                coefficient = -coefficient
            if coefficient:
                for u, v in edges:
                    add(i * n_sites + u, j * n_sites + v, coefficient)
                    add(i * n_sites + v, j * n_sites + u, coefficient)
    for i in range(n_res):
        for u in range(n_sites):
            add(i * n_sites + u, i * n_sites + u, -l1)
            for v in range(u + 1, n_sites):
                add(i * n_sites + u, i * n_sites + v, 2 * l1)
    for u in range(n_sites):
        for i in range(n_res):
            for j in range(i + 1, n_res):
                add(i * n_sites + u, j * n_sites + u, l2)
    adjacent = {edge for edge in edges}
    for i in range(n_res - 1):
        for u in range(n_sites):
            for v in range(n_sites):
                if u != v and (min(u, v), max(u, v)) not in adjacent:
                    add(i * n_sites + u, (i + 1) * n_sites + v, l3)

    metadata = {
        "sequence": sequence, "penalties": [l1, l2, l3],
        "contact_model": model, "contact_alphabet": alphabet,
        "contact_values": contact["values"],
        "lattice_dimensions": list(dims), "site_order": order,
        "coordinates": coords, "lattice_edges": edges,
        "variable_order": "residue-major: index = residue * sites + site",
    }
    return {key: value for key, value in q.items() if value != 0}, n_res * n_sites, l1 * n_res, metadata


def energy(coefficients, offset, bits):
    if any(bit not in "01" for bit in bits):
        raise ValueError("bits must be binary")
    return offset + sum(value * (bits[i] == "1") * (bits[j] == "1")
                        for (i, j), value in coefficients.items())


def validate(bits, metadata):
    """Independently evaluate fold validity and signed contact energy."""
    sequence = metadata["sequence"]
    sites = len(metadata["coordinates"])
    if len(bits) != len(sequence) * sites or any(bit not in "01" for bit in bits):
        raise ValueError("bitstring length or alphabet is invalid")
    placements = [[u for u in range(sites) if bits[i * sites + u] == "1"]
                  for i in range(len(sequence))]
    if any(len(row) != 1 for row in placements):
        return False, None, None
    path = [row[0] for row in placements]
    edges = {tuple(edge) for edge in metadata["lattice_edges"]}
    adjacent = lambda a, b: (min(a, b), max(a, b)) in edges
    if len(set(path)) != len(path) or any(not adjacent(a, b) for a, b in zip(path, path[1:])):
        return False, None, None
    indices = {aa: i for i, aa in enumerate(metadata["contact_alphabet"])}
    values = metadata["contact_values"]
    contact_energy = sum(values[indices[sequence[i]]][indices[sequence[j]]]
                         for i in range(len(path)) for j in range(i + 2, len(path))
                         if adjacent(path[i], path[j]))
    if metadata["contact_model"] == "toy":
        contact_energy = -contact_energy
    return True, contact_energy, path
