"""Published 48-residue instances and an independent coordinate scorer."""

from .model import load_mj85_table_vi


PAPER_URL = "https://juser.fz-juelich.de/record/1047045/files/8n7p-7lh2.pdf"
ORIGINAL_SEQUENCES_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC2242573/"
DIMENSIONS = (4, 4, 3)
REFERENCE_STATUS = "reported_exhaustively_established"

BENCHMARKS = {
    "irback2025_mj1": ("A", "FRTRPLNHDFYNYKIWEPFKPADFPKAWDRMLDHVWDSMASWGHQHCS", -2585),
    "irback2025_mj2": ("A", "CDLPPFTYRHHGNDFWKNYEMIKHWDLWRDMFRAFWSDPVKASPHQAS", -2592),
    "irback2025_mj3": ("A", "FRTPWVSHQFYAYKLMEHFKWGDFCRNMDKWIDSLPDRWNPAPHDHAS", -2609),
    "irback2025_mj4": ("B", "KDKIHFRMNYGYPAWDAQSVKDLTCPRDWHFPHMRDPSHNWELAFFWS", -2587),
    "irback2025_mj5": ("B", "ENDVTMDMDPSPCLFRIHNLPRAHSFDRFGWHQFDKYHYKWKWAWAPS", -2615),
    "irback2025_mj6": ("B", "EHDAQLDFDWSRWTWHGRNSYHAPAMYRWPVHDMDKPNPKFKIFFLCS", -2624),
}


def definition(identifier):
    topology, sequence, minimum_hundredths = BENCHMARKS[identifier]
    if len(sequence) != 48 or any(aa not in "ACDEFGHIKLMNPQRSTVWY" for aa in sequence):
        raise ValueError("benchmark sequence is not 48 standard residues")
    return {
        "id": identifier, "topology": topology, "sequence": sequence,
        "lattice_dimensions": DIMENSIONS, "site_order": "x-fastest",
        "contact_model": "MJ85_TableVI", "full_occupancy": True,
        "reference": {
            "category": "feasible_fold_minimum",
            "status": REFERENCE_STATUS,
            "energy": minimum_hundredths / 100,
            "source": PAPER_URL,
            "source_location": "Table I; exhaustive enumeration in Section II F",
            "ground_state_coordinates": None,
        },
    }


def evaluate_coordinates(identifier, coordinates):
    """Score a complete Hamiltonian path without using QUBO metadata or edges."""
    benchmark = definition(identifier)
    positions = tuple(tuple(point) for point in coordinates)
    if len(positions) != 48 or any(len(point) != 3 for point in positions):
        raise ValueError("exactly 48 three-dimensional coordinates are required")
    sites = {(x, y, z) for x in range(4) for y in range(4) for z in range(3)}
    if any(any(type(v) is not int for v in point) for point in positions):
        raise ValueError("coordinates must be integers")
    if set(positions) != sites or len(set(positions)) != 48:
        raise ValueError("each in-bounds lattice site must be occupied exactly once")
    if any(sum(abs(a - b) for a, b in zip(p, q)) != 1
           for p, q in zip(positions, positions[1:])):
        raise ValueError("consecutive residues must be face-adjacent")
    position_to_residue = {point: i for i, point in enumerate(positions)}
    matrix = load_mj85_table_vi()
    indices = {aa: i for i, aa in enumerate(matrix["alphabet"])}
    sequence = benchmark["sequence"]
    contacts = []
    hundredths = 0
    for x, y, z in sites:
        for neighbor in ((x + 1, y, z), (x, y + 1, z), (x, y, z + 1)):
            if neighbor not in sites:
                continue
            i, j = position_to_residue[x, y, z], position_to_residue[neighbor]
            if abs(i - j) == 1:
                continue
            i, j = sorted((i, j))
            contacts.append((i, j))
            hundredths += round(100 * matrix["values"][indices[sequence[i]]][indices[sequence[j]]])
    if len(contacts) != 57 or len(set(contacts)) != 57:
        raise AssertionError("full-grid fold must have 57 distinct nonbonded contacts")
    return {"contact_count": 57, "contacts": tuple(sorted(contacts)),
            "energy_hundredths": hundredths, "energy": hundredths / 100}


def reference_comparison(identifier, coordinates):
    """A valid fold can match, exceed, or fall below the published value."""
    scored = evaluate_coordinates(identifier, coordinates)
    reference = round(100 * definition(identifier)["reference"]["energy"])
    scored["reference_relation"] = (
        "matches_published_reference_energy" if scored["energy_hundredths"] == reference
        else "below_published_reference_investigate" if scored["energy_hundredths"] < reference
        else "above_published_reference_energy"
    )
    return scored
