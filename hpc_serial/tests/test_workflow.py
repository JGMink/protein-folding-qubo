"""Cheap checks for protected runs and summary grouping; no research budgets."""
import csv
import io
import json
from pathlib import Path
import os
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]

class WorkflowTests(unittest.TestCase):
    def test_run_directory_and_summary(self):
        with tempfile.TemporaryDirectory() as directory:
            environment = dict(os.environ, RESULTS_DIR=directory)
            command = ["bash", str(ROOT / "run/run_experiment.sh"), "fold8", "check",
                       "--sweeps", "2", "--trials", "2", "--seed", "7", "--checkpoint", "1"]
            subprocess.run(command, env=environment, check=True, capture_output=True)
            run = Path(directory) / "check"
            log = run / "trials.jsonl"
            original = log.read_bytes()
            records = [json.loads(line) for line in original.splitlines()]
            self.assertEqual([r["seed"] for r in records], [7, 8])
            self.assertEqual([r["trial"] for r in records], [0, 1])
            self.assertTrue(all(r["attempted_updates"] == 128 for r in records))
            self.assertEqual((run / "status.txt").read_text().strip(), "complete")
            self.assertTrue((run / "checksums.sha256").exists())
            self.assertNotEqual(subprocess.run(command, env=environment, capture_output=True).returncode, 0)
            self.assertEqual(log.read_bytes(), original)
            # Equal temperatures/budgets must not merge different update modes.
            other = Path(directory) / "other.jsonl"
            records[0]["selection"] = "random"
            other.write_text(json.dumps(records[0]) + "\n")
            summary = subprocess.check_output([sys.executable, str(ROOT / "analysis/summarize.py"),
                                               str(log), str(other)], text=True)
            rows = list(csv.DictReader(io.StringIO(summary)))
            self.assertEqual(len(rows), 2)
            self.assertEqual({row["selection"] for row in rows}, {"random", "sequential"})
            self.assertEqual(sorted(int(row["trials"]) for row in rows), [1, 2])
            self.assertTrue(all("mean_final_energy" in row for row in rows))

if __name__ == "__main__":
    unittest.main()
