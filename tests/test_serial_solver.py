"""Portable QUBO format, notebook equivalence, and compiled-solver checks."""

import ast
from collections import defaultdict
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from qubo_bench.export import PRESETS, export
from qubo_bench.io import read_qubo, write_qubo
from qubo_bench.model import DATA, build, energy, validate


def notebook_builder(cell):
    cells = json.loads((ROOT / "notebooks/protein_folding_experiments.ipynb").read_text(encoding="utf-8"))["cells"]
    tree = ast.parse("".join(cells[cell]["source"]))
    functions = ast.Module(body=[node for node in tree.body if isinstance(node, ast.FunctionDef)], type_ignores=[])
    scope = {"np": np, "defaultdict": defaultdict, "MJ_MATRIX": np.array(DATA["signed20"]["values"]),
             "AA_TO_IDX": {aa: i for i, aa in enumerate(DATA["signed20"]["alphabet"] )}}
    exec(compile(functions, "<notebook>", "exec"), scope)
    return scope["build_protein_qubo"]


@pytest.fixture(scope="session")
def solver(tmp_path_factory):
    supplied = os.environ.get("QUBO_SA_EXE")
    if supplied:
        return Path(supplied).resolve()
    compiler = shutil.which("c++") or shutil.which("g++") or shutil.which("clang++")
    if not compiler:
        pytest.skip("set QUBO_SA_EXE or install a C++17 compiler")
    path = tmp_path_factory.mktemp("solver") / "qubo_sa"
    subprocess.run([compiler, "-std=c++17", "-O2", "-Wall", "-Wextra", "-o", str(path),
                    str(ROOT / "solver/qubo_sa.cpp")], check=True)
    return path


def call(solver, path, *args):
    run = subprocess.run([str(solver), str(path), *map(str, args)], check=True, capture_output=True, text=True)
    return [json.loads(line) for line in run.stdout.splitlines()]


def test_export_roundtrip_and_witnesses(tmp_path):
    for name, (sequence, dims, order, model, penalties, category, value, witness) in PRESETS.items():
        path = tmp_path / f"{name}.qubo"
        export(path, sequence, dims, order, model, penalties, (category, value, witness))
        coefficients, n, offset, metadata = read_qubo(path)
        assert n == len(sequence) * np.prod(dims)
        assert metadata["variable_order"].startswith("residue-major")
        assert metadata["penalties"] == list(map(float, penalties))
        assert metadata["reference"]["category"] == category
        if witness:
            bits = ["0"] * n
            for i, site in enumerate(witness):
                bits[i * len(metadata["coordinates"]) + site] = "1"
            valid, direct, path_result = validate("".join(bits), metadata)
            assert valid and path_result == list(witness)
            assert direct == pytest.approx(value)
            assert energy(coefficients, offset, "".join(bits)) == pytest.approx(value)
        again = tmp_path / f"again-{name}.qubo"
        write_qubo(again, coefficients, n, offset, metadata)
        assert again.read_bytes() == path.read_bytes()
        assert path.read_bytes() == (ROOT / "benchmarks" / f"{name}.qubo").read_bytes()


def test_extracted_contact_table_matches_notebooks():
    for name, cell in (("protein_folding_experiments.ipynb", 7), ("protein_folding_qubo.ipynb", 47)):
        cells = json.loads((ROOT / "notebooks" / name).read_text(encoding="utf-8"))["cells"]
        tree = ast.parse("".join(cells[cell]["source"]))
        target = "MJ_MATRIX" if cell == 7 else "MJ_MATRIX_RAW"
        assignments = ast.Module(body=[node for node in tree.body if isinstance(node, ast.Assign)
                                       and any(isinstance(t, ast.Name) and t.id == target for t in node.targets)],
                                 type_ignores=[])
        scope = {"np": np}
        exec(compile(assignments, "<notebook matrix>", "exec"), scope)
        np.testing.assert_array_equal(DATA["signed20"]["values"], scope[target])


def test_exported_toy_has_proven_complete_qubo_minimum():
    coefficients, n, offset, metadata = read_qubo(ROOT / "benchmarks/toy16.qubo")
    states = np.arange(1 << n, dtype=np.uint32)
    energies = np.full(len(states), offset)
    for (i, j), value in coefficients.items():
        energies += value * ((states >> i) & 1) * ((states >> j) & 1)
    assert energies.min() == -1
    minimizers = np.flatnonzero(energies == -1)
    assert len(minimizers) == 8
    for state in minimizers:
        bits = "".join(str((int(state) >> i) & 1) for i in range(n))
        valid, direct, _ = validate(bits, metadata)
        assert valid and direct == -1


@pytest.mark.parametrize("cell,name", [(17, "fold64"), (34, "fold576"), (57, "witness1024")])
def test_sparse_export_matches_notebook_builder(cell, name):
    sequence, dims, order, model, penalties, *_ = PRESETS[name]
    coefficients, n, offset, metadata = build(sequence, dims, penalties, model, order)
    coords = metadata["coordinates"]
    adj = np.zeros((len(coords), len(coords)), dtype=int)
    for u, v in metadata["lattice_edges"]:
        adj[u, v] = adj[v, u] = 1
    linear, quadratic, notebook_offset = notebook_builder(cell)(sequence, adj, *penalties)
    expected = defaultdict(float)
    for i, value in linear.items():
        expected[i, i] += value
    for key, value in quadratic.items():
        expected[key] += value
    assert offset == pytest.approx(notebook_offset)
    assert set(coefficients) == {key for key, value in expected.items() if value}
    for key, value in coefficients.items():
        assert value == pytest.approx(expected[key], abs=1e-10)


def test_3d_export_matches_tutorial_builder():
    sequence = "CFIK"
    coefficients, _, offset, metadata = build(sequence, (2, 2, 2), (1.5, 2, 2), "signed20", "x-fastest")
    adj = np.zeros((8, 8), dtype=int)
    for u, v in metadata["lattice_edges"]:
        adj[u, v] = adj[v, u] = 1
    cells = json.loads((ROOT / "notebooks/protein_folding_qubo.ipynb").read_text(encoding="utf-8"))["cells"]
    tree = ast.parse("".join(cells[53]["source"]))
    functions = ast.Module(body=[node for node in tree.body if isinstance(node, ast.FunctionDef)], type_ignores=[])
    scope = {"np": np, "defaultdict": defaultdict, "MJ_MATRIX_RAW": np.array(DATA["signed20"]["values"]),
             "AA_TO_IDX": {aa: i for i, aa in enumerate(DATA["signed20"]["alphabet"])}}
    exec(compile(functions, "<tutorial builder>", "exec"), scope)
    linear, quadratic, notebook_offset = scope["build_protein_qubo"](sequence, adj, 1.5, 2, 2)
    expected = defaultdict(float)
    for i, value in linear.items():
        expected[i, i] += value
    for key, value in quadratic.items():
        expected[key] += value
    assert offset == notebook_offset
    assert set(coefficients) == {key for key, value in expected.items() if value}
    for key, value in coefficients.items():
        assert value == pytest.approx(expected[key])


def test_reject_malformed_format(tmp_path):
    path = tmp_path / "bad.qubo"
    path.write_text('QUBO1\nn 2\noffset 0\nmeta {}\nq 0 1 1\nq 0 1 2\n')
    with pytest.raises(ValueError, match="duplicate"):
        read_qubo(path)


def test_cpp_energy_and_flip_delta(solver):
    rng = random.Random(11)
    for name in ("toy16", "fold64", "fold576", "witness1024"):
        path = ROOT / "benchmarks" / f"{name}.qubo"
        coefficients, n, offset, _ = read_qubo(path)
        for _ in range(3):
            bits = "".join(rng.choices("01", k=n))
            i = rng.randrange(n)
            flipped = bits[:i] + ("0" if bits[i] == "1" else "1") + bits[i + 1:]
            direct = energy(coefficients, offset, bits)
            delta = energy(coefficients, offset, flipped) - direct
            assert call(solver, path, "--energy", bits)[0]["energy"] == pytest.approx(direct, abs=1e-7)
            assert call(solver, path, "--delta", bits, i)[0]["delta"] == pytest.approx(delta, abs=1e-7)


def test_cpp_metropolis_schedule_counts_and_reproducibility(solver, tmp_path):
    path = tmp_path / "one.qubo"
    write_qubo(path, {(0, 0): 1}, 1, 2.5)
    downhill = call(solver, path, "--initial", "1", "--schedule", "fixed", "--temp-start", 1,
                    "--sweeps", 1)[0]
    assert downhill["best_energy"] == downhill["final_energy"] == 2.5
    assert downhill["accepted_updates"] == 1
    uphill = [call(solver, path, "--initial", "0", "--schedule", "fixed", "--temp-start", 1,
                   "--sweeps", 1, "--seed", seed)[0]["accepted_updates"] for seed in range(200)]
    assert 0.25 < sum(uphill) / len(uphill) < 0.5  # exp(-1) ~= 0.368
    result = call(solver, path, "--initial", "0", "--schedule", "linear", "--temp-start", 2,
                  "--temp-end", 0.5, "--sweeps", 5, "--checkpoint", 2, "--trials", 2,
                  "--seed", 7, "--target", 2.5)
    repeat = call(solver, path, "--initial", "0", "--schedule", "linear", "--temp-start", 2,
                  "--temp-end", 0.5, "--sweeps", 5, "--checkpoint", 2, "--trials", 2,
                  "--seed", 7, "--target", 2.5)
    for first, second in zip(result, repeat):
        for key in ("best_bits", "best_energy", "final_bits", "final_energy", "accepted_updates"):
            assert first[key] == second[key]
        assert first["attempted_updates"] == 5
        assert [p["sweep"] for p in first["progress"]] == [0, 2, 4, 5]
        assert first["first_hit"]["sweep"] == 0
    zero = call(solver, path, "--initial", "1", "--sweeps", 0)[0]
    assert zero["attempted_updates"] == 0 and zero["best_energy"] == 3.5
    assert len(zero["progress"]) == 1


def test_cpp_exact_tiny_and_large_execution(solver, tmp_path):
    path = tmp_path / "tiny.qubo"
    write_qubo(path, {(0, 0): -1, (1, 1): -2, (0, 1): 0.5}, 2, 3)
    run = call(solver, path, "--initial", "00", "--schedule", "fixed", "--temp-start", 0.001,
               "--sweeps", 2)[0]
    assert run["best_bits"] == "11" and run["best_energy"] == pytest.approx(0.5)
    large = ROOT / "benchmarks/witness1024.qubo"
    result = call(solver, large, "--sweeps", 2, "--seed", 4)[0]
    assert result["attempted_updates"] == 2048
    assert len(result["best_bits"]) == 1024
    coefficients, _, offset, _ = read_qubo(large)
    assert result["best_energy"] == pytest.approx(energy(coefficients, offset, result["best_bits"]), abs=1e-6)
