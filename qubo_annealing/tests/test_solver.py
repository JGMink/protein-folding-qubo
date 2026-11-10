"""Standalone checks of fixed inputs, solver arithmetic, and saved folds."""

import csv
import hashlib
import io
import json
from pathlib import Path
import random
import subprocess
import sys
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
        cls.executable = ROOT / "qubo_anneal"

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

    # Two residues on two adjacent sites; variable = residue * 2 + site, so 1001 and 0110 are the
    # only folds (energy 0). `meta` is the fourth line of the file.
    TWO_RESIDUE_BODY = ('QUBO1\nn 4\noffset 2\nMETA\n'
                        'q 0 0 -1\nq 1 1 -1\nq 2 2 -1\nq 3 3 -1\n'
                        'q 0 1 2\nq 2 3 2\nq 0 2 1\nq 1 3 1\n')

    def two_residue_file(self, directory, name, reference=None):
        meta = {"sequence": "AB", "lattice_edges": [[0, 1]]}
        if reference is not None:
            meta["reference"] = reference
        path = Path(directory) / name
        path.write_text(self.TWO_RESIDUE_BODY.replace("META", "meta " + json.dumps(meta)), encoding="utf-8")
        return path

    def test_checkpoints_option_spaces_samples_evenly(self):
        path = ROOT / "benchmarks/fold8.qubo"
        for schedule in ("fixed", "linear", "geometric"):
            record = self.call(path, "--schedule", schedule, "--temp-start", 1, "--temp-end", 1,
                               "--sweeps", 10, "--checkpoints", 5)[0]
            self.assertEqual(record["checkpoint"], 2)
            self.assertEqual([point["sweep"] for point in record["progress"]], [0, 2, 4, 6, 8, 10])
        for bad in (("--sweeps", 10, "--checkpoints", 3),                  # does not divide evenly
                    ("--sweeps", 10, "--checkpoints", 0),
                    ("--sweeps", 10, "--checkpoints", 5, "--checkpoint", 2),
                    ("--schedule", "paper25", "--checkpoints", 5)):
            with self.assertRaises(subprocess.CalledProcessError):
                self.call(path, *bad)

    def test_fold_metrics_deterministic_two_residue_case(self):
        # A near-zero temperature makes every move deterministic. From 0000 the state after
        # sweep 1 is 1010 (both residues on site 0), after sweep 2 it is 0110, the first fold
        # (reached on flip attempt 6), and sweep 3 changes nothing.
        with tempfile.TemporaryDirectory() as directory:
            path = self.two_residue_file(directory, "two.qubo")
            options = ("--initial", "0000", "--schedule", "fixed", "--temp-start", 0.001,
                       "--sweeps", 3, "--checkpoints", 3)
            record = self.call(path, *options, "--fold-metrics")[0]
            plain = self.call(path, *options)[0]
        self.assertNotIn("valid_history", plain)
        self.assertNotIn("valid", plain["progress"][0])
        history = record["valid_history"]
        self.assertEqual(len(history), 1)
        self.assertEqual((history[0]["energy"], history[0]["sweep"], history[0]["attempted_updates"]), (0, 2, 6))
        self.assertEqual(record["first_valid"], history[0])
        self.assertEqual(record["first_best_valid"], history[0])
        self.assertEqual(record["best_valid_energy"], 0)
        self.assertEqual(record["best_valid_bits"], "0110")
        self.assertTrue(record["best_is_valid"] and record["final_is_valid"])
        points = record["progress"]
        self.assertEqual([p["sweep"] for p in points], [0, 1, 2, 3])
        self.assertEqual([p["current_energy"] for p in points], [2, 1, 0, 0])
        self.assertEqual([p["valid"] for p in points], [False, False, True, True])
        self.assertEqual([p["best_valid_energy"] for p in points], [None, None, 0, 0])
        self.assertNotIn("acceptance_rate", points[0])
        self.assertEqual([p["acceptance_rate"] for p in points[1:]], [0.5, 0.5, 0.0])
        self.assertEqual([p["mean_energy"] for p in points[1:]], [1, 0, 0])
        self.assertEqual([p["valid_fraction"] for p in points[1:]], [0, 1, 1])
        # Measuring must not change the trajectory.
        for field in ("attempted_updates", "accepted_updates", "best_bits", "final_bits",
                      "best_energy", "final_energy"):
            self.assertEqual(record[field], plain[field])

    def test_fold_metrics_start_states_and_missing_metadata(self):
        path = ROOT / "benchmarks/fold8.qubo"
        n, _, metadata, _ = read_qubo(path)
        sites = len(metadata["coordinates"])
        reference = ["0"] * n
        for residue, site in enumerate(metadata["reference"]["path"]):
            reference[residue * sites + site] = "1"
        reference = "".join(reference)
        # A known valid fold as the starting state is recorded at sweep 0.
        record = self.call(path, "--initial", reference, "--sweeps", 0, "--fold-metrics")[0]
        self.assertEqual(len(record["valid_history"]), 1)
        self.assertEqual(record["valid_history"][0]["sweep"], 0)
        self.assertAlmostEqual(record["best_valid_energy"], -2.84, places=6)
        self.assertEqual(record["best_valid_bits"], reference)
        self.assertTrue(record["progress"][0]["valid"])
        # All zeros is not a fold, so nothing valid has been seen.
        record = self.call(path, "--initial", "0" * n, "--sweeps", 0, "--fold-metrics")[0]
        self.assertEqual(record["valid_history"], [])
        self.assertIsNone(record["best_valid_energy"])
        self.assertIsNone(record["first_valid"])
        self.assertFalse(record["best_is_valid"])
        self.assertFalse(record["progress"][0]["valid"])
        # Without lattice metadata the request fails with an error instead of guessing.
        with tempfile.TemporaryDirectory() as directory:
            bare = Path(directory) / "bare.qubo"
            bare.write_text("QUBO1\nn 1\noffset 0\nmeta {}\nq 0 0 1\n", encoding="utf-8")
            with self.assertRaises(subprocess.CalledProcessError):
                self.call(bare, "--fold-metrics")

    def test_fold_metrics_match_independent_python_validation(self):
        path = ROOT / "benchmarks/fold8.qubo"
        n, offset, metadata, terms = read_qubo(path)
        proven_minimum = metadata["reference"]["energy"]
        found = 0
        for temperature_end in (0.1, 3.0):  # a cold run usually folds; a hotter one often does not
            for record in self.call(path, "--sweeps", 300, "--trials", 3, "--seed", 5, "--temp-start", 10,
                                    "--temp-end", temperature_end, "--checkpoints", 6, "--fold-metrics"):
                self.assertEqual([p["sweep"] for p in record["progress"]], [0, 50, 100, 150, 200, 250, 300])
                for prefix in ("best", "final"):
                    self.assertEqual(record[f"{prefix}_is_valid"],
                                     validate_fold(record[f"{prefix}_bits"], metadata)[0])
                history = record["valid_history"]
                energies = [entry["energy"] for entry in history]
                self.assertEqual(energies, sorted(set(energies), reverse=True))  # strictly decreasing
                # The tracker must not miss a valid fold that the final states prove exists.
                if record["best_is_valid"]:
                    self.assertAlmostEqual(record["best_valid_energy"], record["best_energy"], places=6)
                if record["final_is_valid"]:
                    self.assertLessEqual(record["best_valid_energy"], record["final_energy"] + 1e-6)
                if not history:
                    self.assertIsNone(record["best_valid_energy"])
                    continue
                found += 1
                self.assertEqual(record["best_valid_energy"], energies[-1])
                valid, contact = validate_fold(record["best_valid_bits"], metadata)
                self.assertTrue(valid)
                self.assertAlmostEqual(contact, record["best_valid_energy"], places=6)
                self.assertAlmostEqual(energy(offset, terms, record["best_valid_bits"]),
                                       record["best_valid_energy"], places=6)
                self.assertGreaterEqual(record["best_valid_energy"], proven_minimum - 1e-6)
        self.assertGreater(found, 0)

    def test_ground_state_analysis_uses_only_proven_references(self):
        script = ROOT / "analysis/ground_state.py"
        references = {"proven": {"category": "feasible_fold_minimum", "energy": 0.0},
                      "witness": {"category": "known_valid_configuration", "energy": 0.0}}
        with tempfile.TemporaryDirectory() as directory:
            rows = {}
            for name, reference in references.items():
                path = self.two_residue_file(directory, f"{name}.qubo", reference)
                log = Path(directory) / f"{name}.jsonl"
                self.call(path, "--initial", "0000", "--schedule", "fixed", "--temp-start", 0.001,
                          "--sweeps", 3, "--fold-metrics", "--output", log)
                result = subprocess.run([sys.executable, str(script), str(log)],
                                        check=True, capture_output=True, text=True)
                rows[name] = list(csv.DictReader(io.StringIO(result.stdout)))[0]
            plain = Path(directory) / "plain.jsonl"
            self.call(path, "--sweeps", 1, "--output", plain)
            rejected = subprocess.run([sys.executable, str(script), str(plain)], capture_output=True, text=True)
        proven, witness = rows["proven"], rows["witness"]
        self.assertEqual((proven["ground_state_found"], proven["ground_state_sweep"],
                          proven["ground_state_attempted_updates"]), ("True", "2", "6"))
        self.assertEqual(float(proven["ground_state_energy"]), 0)
        self.assertEqual(proven["reference_category"], "feasible_fold_minimum")
        self.assertEqual((proven["valid_found"], proven["first_valid_sweep"]), ("True", "2"))
        self.assertEqual(proven["below_reference_investigate"], "False")
        # A witness is not proven optimal: valid folds are still reported, but no ground-state claim.
        self.assertEqual(witness["valid_found"], "True")
        self.assertEqual(witness["ground_state_energy"], "")
        self.assertEqual(witness["ground_state_found"], "")
        self.assertNotEqual(rejected.returncode, 0)  # a log without fold metrics is refused



if __name__ == "__main__":
    unittest.main()
