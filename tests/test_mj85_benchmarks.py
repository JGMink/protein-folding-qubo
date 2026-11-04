"""Independent input, geometry, and scoring gates for the published MJ-48 set."""

import hashlib
import json
from pathlib import Path

import pytest

from qubo_bench.experiment import annotate
from qubo_bench.export import export
from qubo_bench.io import read_qubo
from qubo_bench.mj85_benchmarks import BENCHMARKS, definition, evaluate_coordinates
from qubo_bench.model import DATA, build, energy, lattice, load_mj85_table_vi, validate


EXPECTED = {
    "irback2025_mj1": ("A", "FRTRPLNHDFYNYKIWEPFKPADFPKAWDRMLDHVWDSMASWGHQHCS", -2585),
    "irback2025_mj2": ("A", "CDLPPFTYRHHGNDFWKNYEMIKHWDLWRDMFRAFWSDPVKASPHQAS", -2592),
    "irback2025_mj3": ("A", "FRTPWVSHQFYAYKLMEHFKWGDFCRNMDKWIDSLPDRWNPAPHDHAS", -2609),
    "irback2025_mj4": ("B", "KDKIHFRMNYGYPAWDAQSVKDLTCPRDWHFPHMRDPSHNWELAFFWS", -2587),
    "irback2025_mj5": ("B", "ENDVTMDMDPSPCLFRIHNLPRAHSFDRFGWHQFDKYHYKWKWAWAPS", -2615),
    "irback2025_mj6": ("B", "EHDAQLDFDWSRWTWHGRNSYHAPAMYRWPVHDMDKPNPKFKIFFLCS", -2624),
}


def snake():
    plane = [(x, y) for y in range(4)
             for x in (range(4) if y % 2 == 0 else range(3, -1, -1))]
    return tuple((x, y, z) for z in range(3)
                 for x, y in (plane if z % 2 == 0 else reversed(plane)))


def test_matrix_full_audit_and_default_unchanged():
    path = Path(__file__).resolve().parents[1] / "qubo_bench" / "mj85_table_vi.json"
    resource = json.loads(path.read_text(encoding="utf-8"))
    expected_hash = path.with_suffix(".sha256").read_text(encoding="ascii").strip()
    assert hashlib.sha256(path.read_bytes()).hexdigest() == expected_hash
    assert resource["alphabet"] == "CMFILVWYAGTSQNEDHRKP"  # Trp precedes Tyr in print.
    assert resource["residue_labels"] == [
        "CYS", "MET", "PHE", "ILE", "LEU", "VAL", "TRP", "TYR", "ALA", "GLY",
        "THR", "SER", "GLN", "ASN", "GLU", "ASP", "HIS", "ARG", "LYS", "PRO",
    ]
    triangle = resource["lower_triangle_hundredths"]
    assert len(triangle) == 20 and sum(map(len, triangle)) == 210
    assert all(len(row) == i + 1 and all(type(x) is int for x in row)
               for i, row in enumerate(triangle))
    values = load_mj85_table_vi()["values"]
    assert len(values) == 20 and all(len(row) == 20 for row in values)
    assert all(values[i][j] == values[j][i] == triangle[max(i, j)][min(i, j)] / 100
               for i in range(20) for j in range(20))
    assert values[0][0] == -1.06 and values[-1][-1] == 0.26
    assert DATA["signed20"]["alphabet"] == "ARNDCQEGHILKMFPSTWYV"
    assert DATA["signed20"]["values"][0][0] == -0.2


def test_sequences_and_reported_references():
    assert BENCHMARKS == EXPECTED
    for identifier, (topology, sequence, hundredths) in EXPECTED.items():
        item = definition(identifier)
        assert len(sequence) == 48 and set(sequence) <= set("ACDEFGHIKLMNPQRSTVWY")
        assert item["topology"] == topology
        assert item["reference"]["energy"] == hundredths / 100
        assert item["reference"]["status"] == "reported_exhaustively_established"
        assert item["reference"]["ground_state_coordinates"] is None


@pytest.mark.parametrize("identifier", list(EXPECTED))
def test_all_six_presets_export_with_mj85(tmp_path, identifier):
    item = definition(identifier)
    target = tmp_path / f"{identifier}.qubo"
    n, terms = export(target, item["sequence"], item["lattice_dimensions"],
                      item["site_order"], item["contact_model"], (3, 4, 4),
                      benchmark=item)
    coefficients, parsed_n, offset, metadata = read_qubo(target)
    assert parsed_n == n == 2304 and len(coefficients) == terms
    assert offset == 144 and metadata["benchmark_id"] == identifier
    assert metadata["contact_model"] == "MJ85_TableVI"
    assert metadata["contact_resource_sha256"] == load_mj85_table_vi()["resource_sha256"]
    assert metadata["reference"]["status"] == "reported_exhaustively_established"
    assert metadata["reference"]["energy"] == EXPECTED[identifier][2] / 100
    assert metadata["penalties"] == [3, 4, 4]
    assert len(metadata["lattice_edges"]) == 104
    assert build("AC", (2, 2))[3]["contact_model"] == "signed20"


def test_full_open_grid_edges_and_boundaries():
    coordinates, edges = lattice((4, 4, 3), "x-fastest")
    assert len(coordinates) == 48 and len(set(coordinates)) == 48
    assert len(edges) == 104 and len(set(edges)) == 104
    assert all(sum(abs(a - b) for a, b in zip(coordinates[u], coordinates[v])) == 1
               for u, v in edges)
    neighbors = {site: set() for site in coordinates}
    for u, v in edges:
        neighbors[coordinates[u]].add(coordinates[v])
        neighbors[coordinates[v]].add(coordinates[u])
    assert neighbors[0, 0, 0] == {(1, 0, 0), (0, 1, 0), (0, 0, 1)}
    assert len(neighbors[1, 1, 1]) == 6
    assert (3, 0, 0) not in neighbors[0, 0, 0]  # No periodic wrapping.


@pytest.mark.parametrize("identifier", list(EXPECTED))
def test_coordinate_scorer_has_57_unique_contacts(identifier):
    for positions in (snake(), tuple(reversed(snake())),
                      tuple((3 - x, y, z) for x, y, z in snake())):
        scored = evaluate_coordinates(identifier, positions)
        assert scored["contact_count"] == 57
        assert len(scored["contacts"]) == len(set(scored["contacts"])) == 57
        assert scored["energy_hundredths"] == round(100 * scored["energy"])


def test_invalid_coordinate_configurations_rejected():
    identifier = "irback2025_mj5"
    path = list(snake())
    variants = (
        path[:-1],
        path[:-1] + [path[0]],  # Overlap and vacancy.
        path[:-1] + [(4, 0, 0)],  # Out of bounds.
        path[:-1] + [(-1, 0, 0)],
        path[:1] + [path[2], path[1]] + path[3:],  # Broken chain.
    )
    for positions in variants:
        with pytest.raises(ValueError):
            evaluate_coordinates(identifier, positions)


def test_mj5_qubo_assignment_matches_independent_score_and_zero_penalties():
    item = definition("irback2025_mj5")
    coefficients, n, offset, metadata = build(item["sequence"], item["lattice_dimensions"],
                                               model=item["contact_model"], order=item["site_order"])
    assert n == 2304 and offset == 144
    positions = snake()
    indices = {tuple(point): i for i, point in enumerate(metadata["coordinates"])}
    bits = ["0"] * n
    for residue, point in enumerate(positions):
        bits[residue * 48 + indices[point]] = "1"
    bits = "".join(bits)
    scored = evaluate_coordinates(item["id"], positions)
    valid, contact, path = validate(bits, metadata)
    assert valid and len(set(path)) == 48
    rows = [[int(bit) for bit in bits[i * 48:(i + 1) * 48]] for i in range(48)]
    e1 = sum((sum(row) - 1) ** 2 for row in rows)
    e2 = sum(sum(rows[i][site] * rows[j][site] for i in range(48)
                 for j in range(i + 1, 48)) for site in range(48))
    edges = {tuple(edge) for edge in metadata["lattice_edges"]}
    e3 = sum((min(a, b), max(a, b)) not in edges for a, b in zip(path, path[1:]))
    assert (e1, e2, e3) == (0, 0, 0)
    assert round(100 * contact) == scored["energy_hundredths"]
    assert round(100 * energy(coefficients, offset, bits)) == scored["energy_hundredths"]
    assert round(100 * energy(coefficients, 0, bits)) == scored["energy_hundredths"] - 14400
    assert scored["energy_hundredths"] > -2615  # This path is a check, not a ground state.
    occupied = [i * 48 + indices[point] for i, point in enumerate(positions)]
    overlap = list(bits)
    overlap[occupied[1]] = "0"
    overlap[1 * 48 + indices[positions[0]]] = "1"
    missing = list(bits)
    missing[occupied[1]] = "0"
    extra = list(bits)
    extra[1 * 48 + indices[positions[0]]] = "1"
    broken = list(bits)
    for residue in (1, 2):
        broken[occupied[residue]] = "0"
        broken[residue * 48 + indices[positions[3 - residue]]] = "1"
    for candidate in (overlap, missing, extra, broken):
        assert validate("".join(candidate), metadata)[0] is False
    metadata["reference"] = item["reference"]
    record = {"solver": "cpp", "best_bits": bits,
              "best_energy": energy(coefficients, offset, bits),
              "final_bits": bits, "final_energy": energy(coefficients, offset, bits)}
    annotated = annotate(record, coefficients, n, offset, metadata)
    assert annotated["matches_published_reference_energy"] is False
    invalid_bits = "0" * n
    invalid_energy = energy(coefficients, offset, invalid_bits)
    metadata["reference"] = dict(item["reference"], energy=invalid_energy + 1)
    invalid_record = {"solver": "cpp", "best_bits": invalid_bits,
                      "best_energy": invalid_energy, "final_bits": invalid_bits,
                      "final_energy": invalid_energy}
    annotated = annotate(invalid_record, coefficients, n, offset, metadata)
    assert annotated["valid"] is False
    assert annotated["matches_published_reference_energy"] is False
