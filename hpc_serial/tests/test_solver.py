"""Standalone checks of fixed inputs, solver arithmetic, and saved folds."""

import hashlib
import json
from pathlib import Path
import random
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS = {
    "fold8": (64, "c044ce6c9690c281a5d792e6c1591630977def90ab06d4f253828a29877d53d8", -2.84,
              "feasible_fold_minimum"),
    "fold24": (576, "614f5b054cb33f2d296d54d9c54f41115fce078a30de9c67aaa623b0cddb0b8c", -10.05,
               "feasible_fold_minimum"),
    "fold32": (1024, "74804084a21b98e021bfe25031073f1c50e3da9a8990c5a1a1b5ed4c3d3d465e", -14.25,
               "known_valid_configuration"),
}


def read_qubo(path):
    """Small copy of the QUBO1 reader used by the original tests."""
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 4 or lines[0] != "QUBO1":
        raise ValueError("expected QUBO1")
    n = int(lines[1].split()[1])
    offset = float(lines[2].split()[1])
    metadata = json.loads(lines[3][5:])
    terms = []
    seen = set()
    for line in lines[4:]:
        label, i, j, coefficient = line.split()
        i, j = int(i), int(j)
        if label != "q" or not 0 <= i <= j < n or (i, j) in seen:
            raise ValueError("invalid QUBO coefficient")
        seen.add((i, j))
        terms.append((i, j, float(coefficient)))
    return n, offset, metadata, terms


def energy(offset, terms, bits):
    return offset + sum(value for i, j, value in terms if bits[i] == "1" and bits[j] == "1")


def validate_fold(bits, metadata):
    """Adapted from qubo_bench.model.validate; uses only embedded metadata."""
    sequence = metadata["sequence"]
    sites = len(metadata["coordinates"])
    if len(bits) != len(sequence) * sites or set(bits) - {"0", "1"}:
        raise ValueError("invalid bitstring")
    placements = [[site for site in range(sites) if bits[i * sites + site] == "1"]
                  for i in range(len(sequence))]
    if any(len(row) != 1 for row in placements):
        return False, None
    path = [row[0] for row in placements]
    edges = {tuple(edge) for edge in metadata["lattice_edges"]}

    def adjacent(a, b):
        return (min(a, b), max(a, b)) in edges

    if len(set(path)) != len(path) or any(not adjacent(a, b) for a, b in zip(path, path[1:])):
        return False, None
    index = {aa: i for i, aa in enumerate(metadata["contact_alphabet"])}
    values = metadata["contact_values"]
    contact = sum(values[index[sequence[i]]][index[sequence[j]]]
                  for i in range(len(path)) for j in range(i + 2, len(path))
                  if adjacent(path[i], path[j]))
    return True, contact


class SerialSolverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["make"], cwd=ROOT, check=True, capture_output=True, text=True)
        cls.executable = ROOT / "serial"

    def call(self, path, *options):
        result = subprocess.run([str(self.executable), str(path), *map(str, options)],
                                cwd=ROOT, check=True, capture_output=True, text=True)
        return [json.loads(line) for line in result.stdout.splitlines()]

    def test_fixed_inputs_and_reference_witnesses(self):
        for name, (expected_n, digest, reference, category) in BENCHMARKS.items():
            with self.subTest(name=name):
                path = ROOT / "benchmarks" / f"{name}.qubo"
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), digest)
                n, offset, metadata, terms = read_qubo(path)
                self.assertEqual(n, expected_n)
                self.assertEqual(metadata["reference"]["category"], category)
                self.assertEqual(metadata["reference"]["energy"], reference)
                sites = len(metadata["coordinates"])
                bits = ["0"] * n
                for residue, site in enumerate(metadata["reference"]["path"]):
                    bits[residue * sites + site] = "1"
                bits = "".join(bits)
                valid, contact = validate_fold(bits, metadata)
                self.assertTrue(valid)
                self.assertAlmostEqual(contact, reference, places=6)
                self.assertAlmostEqual(energy(offset, terms, bits), reference, places=6)

    def test_cpp_energy_and_incremental_flip_delta(self):
        rng = random.Random(11)
        for name in BENCHMARKS:
            path = ROOT / "benchmarks" / f"{name}.qubo"
            n, offset, _, terms = read_qubo(path)
            bits = "".join(rng.choices("01", k=n))
            index = rng.randrange(n)
            flipped = bits[:index] + ("0" if bits[index] == "1" else "1") + bits[index + 1:]
            direct = energy(offset, terms, bits)
            delta = energy(offset, terms, flipped) - direct
            with self.subTest(name=name):
                self.assertAlmostEqual(self.call(path, "--energy", bits)[0]["energy"], direct, places=6)
                self.assertAlmostEqual(self.call(path, "--delta", bits, index)[0]["delta"], delta, places=6)

    def test_qubo_utilities_link_without_serial_main(self):
        source = r'''
#include "qubo.h"
#include <cmath>
#include <vector>

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    const Qubo q = load_qubo(argv[1]);
    std::vector<unsigned char> first(q.n, 0), second(q.n, 0);
    first[0] = 1;
    second[1] = 1;
    auto first_fields = local_fields(q, first);
    const auto second_fields = local_fields(q, second);
    const double before = energy(q, first), other_before = energy(q, second);
    const double delta = flip_delta(first, first_fields, 2);
    const int change = 1 - 2 * int(first[2]);
    update_fields_after_flip(q, first_fields, 2, change);
    first[2] ^= 1;
    if (std::abs(energy(q, first) - (before + delta)) > 1e-9) return 3;
    const auto recalculated = local_fields(q, first);
    for (int i = 0; i < q.n; ++i)
        if (std::abs(first_fields[i] - recalculated[i]) > 1e-9) return 4;
    if (second_fields != local_fields(q, second) || energy(q, second) != other_before) return 5;
    return 0;
}
'''
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            (path / "utility_check.cpp").write_text(source, encoding="utf-8")
            executable = path / "utility_check"
            subprocess.run(["c++", "-std=c++17", "-O2", "-I", str(ROOT),
                            str(path / "utility_check.cpp"), str(ROOT / "qubo.cpp"),
                            "-o", str(executable)], check=True, capture_output=True, text=True)
            subprocess.run([str(executable), str(ROOT / "benchmarks/fold8.qubo")], check=True)

    def test_execution_log_reproducibility_and_saved_fold_energy(self):
        reproducible = ("seed", "trial", "best_bits", "best_energy", "final_bits",
                        "final_energy", "accepted_updates", "first_best", "target_reached")
        for name, (n, _, _, _) in BENCHMARKS.items():
            path = ROOT / "benchmarks" / f"{name}.qubo"
            _, offset, metadata, terms = read_qubo(path)
            options = ("--sweeps", 2, "--trials", 2, "--seed", 7, "--schedule", "linear",
                       "--temp-start", 2, "--temp-end", 0.5, "--checkpoint", 1)
            first = self.call(path, *options)
            second = self.call(path, *options)
            with self.subTest(name=name):
                self.assertEqual(len(first), 2)
                for a, b in zip(first, second):
                    for field in reproducible:
                        if field == "first_best":
                            self.assertEqual(a[field]["sweep"], b[field]["sweep"])
                            self.assertEqual(a[field]["attempted_updates"], b[field]["attempted_updates"])
                        else:
                            self.assertEqual(a[field], b[field])
                    self.assertEqual(a["benchmark"], name)
                    self.assertEqual(a["attempted_updates"], 2 * n)
                    self.assertEqual(a["sweeps"], 2)
                    self.assertEqual(a["checkpoint"], 1)
                    self.assertIsNone(a["initial_bits"])
                    self.assertIsNone(a["target_energy"])
                    self.assertIsNone(a["first_hit"])
                    self.assertGreaterEqual(a["runtime_seconds"], a["first_best"]["elapsed_seconds"])
                    for prefix in ("best", "final"):
                        bits = a[f"{prefix}_bits"]
                        self.assertEqual(len(bits), n)
                        self.assertAlmostEqual(a[f"{prefix}_energy"], energy(offset, terms, bits), places=6)
                        valid, contact = validate_fold(bits, metadata)
                        if valid:
                            self.assertAlmostEqual(contact, a[f"{prefix}_energy"], places=6)

    def test_first_best_and_optional_target(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "one.qubo"
            path.write_text('QUBO1\nn 1\noffset 2.5\nmeta {}\nq 0 0 1\n', encoding="utf-8")
            record = self.call(path, "--initial", "1", "--schedule", "fixed",
                               "--temp-start", 1, "--sweeps", 1, "--target", 2.5)[0]
            self.assertEqual(record["initial_bits"], "1")
            self.assertEqual(record["target_energy"], 2.5)
            self.assertTrue(record["target_reached"])
            self.assertEqual(record["best_energy"], 2.5)
            self.assertEqual(record["first_best"]["sweep"], 1)
            self.assertEqual(record["first_best"]["attempted_updates"], 1)
            self.assertGreaterEqual(record["first_best"]["elapsed_seconds"], 0)
            self.assertEqual(record["first_hit"]["attempted_updates"], 1)

    def test_paper_stages_random_selection_and_checkpoint_energy(self):
        path = ROOT / 'benchmarks/fold8.qubo'
        n, offset, _, terms = read_qubo(path)
        options = ('--schedule', 'paper25', '--selection', 'random',
                   '--sweeps-per-temperature', 3, '--seed', 19)
        a = self.call(path, *options)[0]
        b = self.call(path, *options)[0]
        self.assertEqual(a['sweeps'], 75)
        self.assertEqual(a['attempted_updates'], 75*n)
        self.assertEqual(a['best_bits'], b['best_bits'])
        self.assertEqual(a['final_bits'], b['final_bits'])
        self.assertEqual(a['accepted_updates'], b['accepted_updates'])
        self.assertEqual(len(a['progress']), 26)
        self.assertAlmostEqual(a['temperature_start'], 1/1.05)
        self.assertAlmostEqual(a['temperature_end'], 1/1.05**25)
        for stage,point in enumerate(a['progress'][1:], 1):
            self.assertEqual(point['sweep'], stage*3)
            self.assertEqual(point['attempted_updates'], stage*3*n)
            self.assertAlmostEqual(point['beta'], 1.05**stage)
            self.assertAlmostEqual(point['current_energy'], energy(offset,terms,point['bits']))
        # An empty polynomial accepts every proposal: with-replacement selection
        # can flip a bit repeatedly, while sequential selection flips every bit.
        with tempfile.TemporaryDirectory() as directory:
            empty=Path(directory)/'empty.qubo'
            empty.write_text('QUBO1\nn 16\noffset 0\nmeta {}\n')
            a=self.call(empty,'--initial','0'*16,'--sweeps',1,'--seed',19,
                        '--selection','random')[0]
            b=self.call(empty,'--initial','0'*16,'--sweeps',1)[0]
            self.assertEqual(a['accepted_updates'],16)
            self.assertNotEqual(a['final_bits'],'1'*16)
            self.assertEqual(b['final_bits'],'1'*16)


if __name__ == "__main__":
    unittest.main()
