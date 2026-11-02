"""Independent small-state checks of the lattice protein QUBO."""

import ast
from collections import defaultdict
import itertools
import json
from pathlib import Path
import sys
import time

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "problem_formation_and_evaluation" / "QUBO_construction"))
sys.path.insert(0, str(ROOT / "problem_formation_and_evaluation" / "energy_calc"))
from qubo_generation import build_E_MJ, build_E1, build_E2, build_E3
from calc_mods import total_energy, is_valid_conformation


def lattice(rows, cols):
    adj = np.zeros((rows * cols, rows * cols), dtype=int)
    for r in range(rows):
        for c in range(cols):
            v = r * cols + c
            for dr, dc in ((1, 0), (0, 1)):
                rr, cc = r + dr, c + dc
                if rr < rows and cc < cols:
                    w = rr * cols + cc
                    adj[v, w] = adj[w, v] = 1
    return adj


def independent_components(chain, b, adj, contacts):
    n, m = b.shape
    contact = -sum(
        contacts.get((chain[i], chain[j]), contacts.get((chain[j], chain[i]), 0))
        for i in range(n) for j in range(i + 2, n)
        for u in range(m) for v in range(m)
        if b[i, u] and b[j, v] and adj[u, v]
    )
    e1 = sum((sum(b[i]) - 1) ** 2 for i in range(n))
    e2 = sum(b[i, u] * b[j, u] for u in range(m)
             for i in range(n) for j in range(i + 1, n))
    e3 = sum(b[i, u] * b[i + 1, v]
             for i in range(n - 1) for u in range(m) for v in range(m)
             if u != v and not adj[u, v])
    return np.array((contact, e1, e2, e3), dtype=float)


def test_exhaustive_binary_states_and_feasibility():
    chain, adj, contacts = "HPCH", lattice(2, 2), {("H", "H"): 1}
    components = (build_E_MJ(chain, adj, contacts), build_E1(chain, 4),
                  build_E2(chain, 4), build_E3(chain, adj))
    rng = np.random.default_rng(7)
    for state in itertools.product((0, 1), repeat=16):
        x = np.array(state)
        b = x.reshape(4, 4)
        expected = independent_components(chain, b, adj, contacts)
        computed = [x @ q @ x + offset for q, _, offset in components]
        np.testing.assert_allclose(computed, expected)
        if rng.random() < 0.002:
            energy, breakdown = total_energy(chain, ''.join(map(str, state)), adj, contacts)
            np.testing.assert_allclose(breakdown, expected)
            assert energy == sum(expected)
        valid = all(sum(b[i]) == 1 for i in range(4)) and all(sum(b[:, u]) <= 1 for u in range(4))
        valid = valid and all(adj[np.argmax(b[i]), np.argmax(b[i + 1])] for i in range(3))
        assert is_valid_conformation(expected) == valid


def test_rectangular_occupancy_and_rejected_malformed_inputs():
    chain, adj = "HH", lattice(1, 3)
    energy, breakdown = total_energy(chain, "100010", adj, {("H", "H"): 1})
    assert energy == 0 and is_valid_conformation(breakdown)
    _, breakdown = total_energy(chain, "100001", adj, {})
    assert breakdown[3] == 1
    for bits in ("10001", "1000100", "10001x"):
        with pytest.raises(ValueError):
            total_energy(chain, bits, adj, {})
    for bad in (np.array([[0, 1], [0, 0]]), np.array([[1, 0], [0, 0]]),
                np.array([[0, 2], [2, 0]])):
        with pytest.raises(ValueError):
            build_E3("HH", bad)
    with pytest.raises(ValueError):
        build_E1('H', 0)
    with pytest.raises(ValueError):
        build_E2('H', 0)
    with pytest.raises(ValueError):
        build_E_MJ('ABBA', lattice(2, 2), {('A', 'B'): 1, ('B', 'A'): 2})


def test_polynomials_match_qubo_coefficients():
    chain, adj = "HPCH", lattice(2, 2)
    for q, poly, offset in (build_E_MJ(chain, adj, {("H", "H"): 1}),
                            build_E1(chain, 4), build_E2(chain, 4), build_E3(chain, adj)):
        reconstructed = np.zeros_like(q)
        for coefficient, i, j in poly:
            reconstructed[min(i, j), max(i, j)] += coefficient
        np.testing.assert_array_equal(reconstructed, q)


def notebook_functions(cell_index):
    cells = json.loads((ROOT / "notebooks" / "protein_folding_experiments.ipynb").read_text(encoding="utf-8"))["cells"]
    tree = ast.parse(''.join(cells[cell_index]['source']))
    functions = ast.Module(body=[node for node in tree.body if isinstance(node, ast.FunctionDef)], type_ignores=[])
    scope = {"np": np, "defaultdict": defaultdict, "MJ_MATRIX": np.array([[-1.0]]),
             "AA_TO_IDX": {"A": 0}}
    exec(compile(functions, '<notebook>', 'exec'), scope)
    return scope


def notebook_matrix(cell_index):
    cells = json.loads((ROOT / "notebooks" / "protein_folding_experiments.ipynb").read_text(encoding="utf-8"))["cells"]
    tree = ast.parse(''.join(cells[cell_index]['source']))
    assignments = ast.Module(body=[node for node in tree.body
                                   if isinstance(node, ast.Assign) and any(
                                       isinstance(target, ast.Name) and target.id == 'MJ_MATRIX'
                                       for target in node.targets)], type_ignores=[])
    scope = {'np': np}
    exec(compile(assignments, '<matrix>', 'exec'), scope)
    return scope['MJ_MATRIX']


@pytest.mark.parametrize("cell", [17, 34, 57])
def test_experiment_qubo_uses_signed_contact_energy(cell):
    scope = notebook_functions(cell)
    adj = lattice(2, 2)
    linear, quadratic, offset = scope['build_protein_qubo']("AAAA", adj, 2, 2, 2)
    path = (0, 1, 3, 2)
    bits = np.zeros(16)
    bits[[i * 4 + u for i, u in enumerate(path)]] = 1
    energy = offset + sum(v * bits[i] for i, v in linear.items())
    energy += sum(v * bits[i] * bits[j] for (i, j), v in quadratic.items())
    assert energy == -1  # Exactly one non-chain contact on a 2x2 Hamiltonian path.


def test_embedded_tutorial_sources_match_scripts():
    cells = json.loads((ROOT / "notebooks" / "protein_folding_qubo.ipynb").read_text(encoding="utf-8"))["cells"]
    for index, relative in ((4, "QUBO_construction/qubo_generation.py"),
                            (5, "energy_calc/calc_mods.py")):
        embedded = ''.join(cells[index]['source']).split('\n', 1)[1]
        actual = (ROOT / "problem_formation_and_evaluation" / relative).read_text(encoding="utf-8")
        assert embedded == actual


def test_experiment_matrices_and_exact_small_benchmarks():
    matrices = [notebook_matrix(index) for index in (7, 28, 50)]
    for mat in matrices:
        assert mat.shape == (20, 20)
        np.testing.assert_array_equal(mat, mat.T)
        np.testing.assert_array_equal(mat, matrices[0])
    alphabet = 'ARNDCQEGHILKMFPSTWYV'
    index = {aa: i for i, aa in enumerate(alphabet)}

    def enumerate_paths(rows, cols, sequence):
        adj = lattice(rows, cols)
        neighbors = [tuple(np.flatnonzero(adj[u])) for u in range(rows * cols)]
        path = []
        count, minimum = 0, float('inf')

        def visit(mask, energy):
            nonlocal count, minimum
            if len(path) == len(sequence):
                count += 1
                minimum = min(minimum, energy)
                return
            for v in neighbors[path[-1]]:
                if mask & (1 << int(v)):
                    continue
                increment = sum(matrices[0][index[sequence[i]], index[sequence[len(path)]]]
                                for i, u in enumerate(path[:-1]) if adj[u, v])
                path.append(int(v))
                visit(mask | (1 << int(v)), energy + increment)
                path.pop()

        for start in range(rows * cols):
            path.append(start)
            visit(1 << start, 0.0)
            path.pop()
        return count, minimum

    count, minimum = enumerate_paths(2, 4, 'CFIKDLWV')
    assert count == 28
    assert minimum == pytest.approx(-2.84)
    count, minimum = enumerate_paths(4, 6, 'ENDVTMDMDPSPCLFRIHNLPRAH')
    assert count == 7220
    assert minimum == pytest.approx(-10.05)
    sequence = 'ENDVTMDMDPSPCLFRIHNLPRAHSFDRFGWH'
    path = [9, 10, 14, 18, 22, 21, 25, 26, 30, 29, 28, 24, 20, 16, 17, 13,
            12, 8, 4, 0, 1, 5, 6, 2, 3, 7, 11, 15, 19, 23, 27, 31]
    adj = lattice(8, 4)
    assert len(set(path)) == 32 and all(adj[u, v] for u, v in zip(path, path[1:]))
    contacts = [(i, j) for i in range(32) for j in range(i + 2, 32)
                if adj[path[i], path[j]]]
    assert len(contacts) == 21
    assert sum(matrices[0][index[sequence[i]], index[sequence[j]]]
               for i, j in contacts) == pytest.approx(-14.25)


def test_full_size_tutorial_builder_and_validator_use_signed_energies():
    cells = json.loads((ROOT / 'notebooks' / 'protein_folding_qubo.ipynb').read_text(encoding='utf-8'))['cells']
    matrix_tree = ast.parse(''.join(cells[47]['source']))
    matrix_code = ast.Module(body=[node for node in matrix_tree.body
                                   if isinstance(node, ast.Assign) and any(
                                       isinstance(target, ast.Name) and target.id == 'MJ_MATRIX_RAW'
                                       for target in node.targets)], type_ignores=[])
    scope = {'np': np, 'defaultdict': defaultdict, 'AA_TO_IDX': {'A': 0}}
    exec(compile(matrix_code, '<tutorial matrix>', 'exec'), scope)
    np.testing.assert_array_equal(scope['MJ_MATRIX_RAW'], notebook_matrix(7))
    for index in (53, 57):
        tree = ast.parse(''.join(cells[index]['source']))
        functions = ast.Module(body=[node for node in tree.body if isinstance(node, ast.FunctionDef)], type_ignores=[])
        exec(compile(functions, '<tutorial>', 'exec'), scope)
    adj = lattice(2, 2)
    linear, quadratic, offset = scope['build_protein_qubo']('AAAA', adj, 2, 2, 2)
    sample = {i: 0 for i in range(16)}
    for i, site in enumerate((0, 1, 3, 2)):
        sample[i * 4 + site] = 1
    valid, energy, violations = scope['validate_solution'](sample, 'AAAA', adj, 4)
    assert valid and all(value == 0 for value in violations.values())
    assert energy == pytest.approx(-0.20)
    qubo = offset + sum(value * sample[i] for i, value in linear.items())
    qubo += sum(value * sample[i] * sample[j] for (i, j), value in quadratic.items())
    assert qubo == pytest.approx(energy)
    with pytest.raises(ValueError):
        scope['validate_solution']({0: 1}, 'AAAA', adj, 4)


@pytest.mark.parametrize('cell,rows,cols,sequence,path,expected', [
    (34, 4, 6, 'ENDVTMDMDPSPCLFRIHNLPRAH',
     [1, 0, 6, 12, 18, 19, 20, 21, 22, 23, 17, 16, 15, 14, 13, 7, 8, 2, 3, 9, 10, 11, 5, 4], -10.05),
    (57, 8, 4, 'ENDVTMDMDPSPCLFRIHNLPRAHSFDRFGWH',
     [9, 10, 14, 18, 22, 21, 25, 26, 30, 29, 28, 24, 20, 16, 17, 13, 12, 8, 4, 0, 1, 5, 6, 2, 3, 7, 11, 15, 19, 23, 27, 31], -14.25),
])
def test_large_experiment_qubo_matches_path_energy(cell, rows, cols, sequence, path, expected):
    scope = notebook_functions(cell)
    scope['MJ_MATRIX'] = notebook_matrix(7)
    scope['AA_TO_IDX'] = {aa: i for i, aa in enumerate('ARNDCQEGHILKMFPSTWYV')}
    if cell == 57:  # This section numbers sites by column, unlike the row-major witness.
        path = [(site % cols) * rows + site // cols for site in path]
        coords = [(r, c) for c in range(cols) for r in range(rows)]
        adj = np.array([[int(abs(r - rr) + abs(c - cc) == 1)
                         for rr, cc in coords] for r, c in coords])
    else:
        adj = lattice(rows, cols)
    linear, quadratic, offset = scope['build_protein_qubo'](sequence, adj, 1.5, 2, 2)
    x = np.zeros(len(sequence) * len(path))
    for i, site in enumerate(path):
        x[i * len(path) + site] = 1
    qubo = offset + sum(value * x[i] for i, value in linear.items())
    qubo += sum(value * x[i] * x[j] for (i, j), value in quadratic.items())
    direct = sum(scope['MJ_MATRIX'][scope['AA_TO_IDX'][sequence[i]], scope['AA_TO_IDX'][sequence[j]]]
                 for i in range(len(path)) for j in range(i + 2, len(path)) if adj[path[i], path[j]])
    assert qubo == pytest.approx(direct)
    assert direct == pytest.approx(expected)


def test_48_residue_tutorial_witness_matches_qubo_and_validator():
    cells = json.loads((ROOT / 'notebooks' / 'protein_folding_qubo.ipynb').read_text(encoding='utf-8'))['cells']
    scope = {'np': np, 'defaultdict': defaultdict}
    for cell in (47, 51, 53, 57):
        tree = ast.parse(''.join(cells[cell]['source']))
        selected = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
        if cell == 47:
            selected = [node for node in tree.body if isinstance(node, ast.Assign)
                        and any(isinstance(target, ast.Name) and target.id in ('AMINO_ACIDS', 'AA_TO_IDX', 'MJ_MATRIX_RAW')
                                for target in node.targets)]
        exec(compile(ast.Module(body=selected, type_ignores=[]), '<tutorial>', 'exec'), scope)
    sequence = 'ENDVTMDMDPSPCLFRIHNLPRAHSFDRFGWHQFDKYHYKWKWAWAPS'
    adj, coords = scope['build_cubic_lattice']((4, 4, 3))
    location = {coord: site for site, coord in enumerate(coords)}
    plane = [(x, y) for y in range(4) for x in (range(4) if y % 2 == 0 else range(3, -1, -1))]
    path = [location[x, y, z] for z in range(3)
            for x, y in (plane if z % 2 == 0 else reversed(plane))]
    assert len(path) == 48 and len(set(path)) == 48
    assert all(adj[u, v] for u, v in zip(path, path[1:]))
    linear, quadratic, offset = scope['build_protein_qubo'](sequence, adj)
    sample = {i: 0 for i in range(48 * 48)}
    for i, site in enumerate(path):
        sample[i * 48 + site] = 1
    valid, direct, penalties = scope['validate_solution'](sample, sequence, adj, 48)
    assert valid and all(value == 0 for value in penalties.values())
    qubo = offset + sum(value * sample[i] for i, value in linear.items())
    qubo += sum(value * sample[i] * sample[j] for (i, j), value in quadratic.items())
    assert qubo == pytest.approx(direct)


def test_neal_samples_match_direct_energy_when_available():
    neal = pytest.importorskip('neal')
    dimod = pytest.importorskip('dimod')
    scope = notebook_functions(17)
    scope['MJ_MATRIX'] = notebook_matrix(7)
    scope['AA_TO_IDX'] = {aa: i for i, aa in enumerate('ARNDCQEGHILKMFPSTWYV')}
    cells = json.loads((ROOT / 'notebooks' / 'protein_folding_experiments.ipynb').read_text(encoding='utf-8'))['cells']
    tree = ast.parse(''.join(cells[19]['source']))
    functions = ast.Module(body=[node for node in tree.body if isinstance(node, ast.FunctionDef)], type_ignores=[])
    exec(compile(functions, '<notebook validator>', 'exec'), scope)
    sequence, adj = 'CFIKDLWV', lattice(2, 4)
    linear, quadratic, offset = scope['build_protein_qubo'](sequence, adj, 3, 4, 4)
    bqm = dimod.BinaryQuadraticModel(linear, quadratic, offset, vartype='BINARY')
    samples = neal.SimulatedAnnealingSampler().sample(
        bqm, num_reads=5, num_sweeps=100, beta_range=(1, 50), seed=42)
    assert len(samples) == 5
    for datum in samples.data():
        valid, contact, penalties, path = scope['validate_solution'](datum.sample, sequence, adj)
        direct = contact + 3 * penalties['E1'] + 4 * penalties['E2'] + 4 * penalties['E3']
        assert datum.energy == pytest.approx(direct, abs=1e-6)
        if valid:
            assert len(set(path)) == len(sequence)


@pytest.mark.parametrize('matrix_cell,func_cells,dimensions,sequence', [
    (28, (30, 32, 34, 36, 40), (4, 6), 'ENDVTMDMDPSPCLFRIHNLPRAH'),
    (50, (52, 55, 57, 59, 61), (8, 4), 'ENDVTMDMDPSPCLFRIHNLPRAHSFDRFGWH'),
])
def test_neal_experiment_processing_when_available(matrix_cell, func_cells, dimensions, sequence):
    neal = pytest.importorskip('neal')
    dimod = pytest.importorskip('dimod')
    cells = json.loads((ROOT / 'notebooks' / 'protein_folding_experiments.ipynb').read_text(encoding='utf-8'))['cells']
    scope = {'np': np, 'defaultdict': defaultdict, 'time': time,
             'BinaryQuadraticModel': dimod.BinaryQuadraticModel,
             'SimulatedAnnealingSampler': neal.SimulatedAnnealingSampler,
             'AA_TO_IDX': {aa: i for i, aa in enumerate('ARNDCQEGHILKMFPSTWYV')},
             'MJ_MATRIX': notebook_matrix(matrix_cell)}
    for cell in func_cells:
        tree = ast.parse(''.join(cells[cell]['source']))
        functions = ast.Module(body=[node for node in tree.body if isinstance(node, ast.FunctionDef)], type_ignores=[])
        exec(compile(functions, '<experiment>', 'exec'), scope)
    built = scope['build_2d_lattice'](*dimensions)
    adj = built[0] if isinstance(built, tuple) else built
    linear, quadratic, offset = scope['build_protein_qubo'](sequence, adj, 3, 4, 4)
    run = scope['run_sa_experiment'](linear, quadratic, offset, 2, 100, (1, 50),
                                     sequence, adj, 3, 4, 4, seed=42)
    assert len(run['all_results']) == 2
    for result in run['all_results']:
        assert result['QUBO_energy'] == pytest.approx(result['E_total'], abs=1e-6)
