"""Summarize serial solver JSONL logs by benchmark and experiment settings."""

import argparse
import csv
import json
from pathlib import Path
from statistics import mean, median
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("logs", nargs="+", type=Path, help="JSONL files, for example results/*.jsonl")
    args = parser.parse_args()
    groups = {}
    for path in args.logs:
        with path.open(encoding="utf-8") as stream:
            for line_number, line in enumerate(stream, 1):
                if not line.strip():
                    continue
                record = json.loads(line)
                try:
                    key = (
                        record["input"], record["benchmark"], record["schedule"],
                        record.get("selection", "sequential"), record.get("sweeps_per_temperature", 0),
                        record["temperature_start"],
                        record["temperature_end"], record["sweeps"], record["checkpoint"],
                        record["initial_bits"], record["target_energy"],
                    )
                    groups.setdefault(key, []).append(record)
                except KeyError as error:
                    parser.error(f"{path}:{line_number}: missing field {error}")
    writer = csv.writer(sys.stdout)
    writer.writerow(("input", "benchmark", "schedule", "selection", "sweeps_per_temperature", "temp_start", "temp_end", "sweeps", "checkpoint",
                     "initial_bits", "target_energy", "trials", "mean_best_energy",
                     "median_best_energy", "best_observed_energy", "mean_final_energy", "mean_runtime_seconds",
                     "target_hits", "target_hit_rate"))
    for key, records in sorted(groups.items(), key=lambda item: tuple(str(v) for v in item[0])):
        energies = [record["best_energy"] for record in records]
        runtimes = [record["runtime_seconds"] for record in records]
        hits = sum(bool(record["target_reached"]) for record in records)
        has_target = key[-1] is not None
        writer.writerow((*key, len(records), f"{mean(energies):.6f}", f"{median(energies):.6f}",
                         f"{min(energies):.6f}",
                         f"{mean(record['final_energy'] for record in records):.6f}", f"{mean(runtimes):.6f}",
                         hits if has_target else "", f"{hits / len(records):.3f}" if has_target else ""))


if __name__ == "__main__":
    main()
