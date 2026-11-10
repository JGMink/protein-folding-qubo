"""Ground-state success per trial, from logs written with --fold-metrics.

The ground state is the proven minimum energy over valid folds. It is read from the
`reference` metadata of each record's input file, and only when that reference is a proven
minimum (feasible_fold_minimum or complete_qubo_global_minimum). Witnesses such as
known_valid_configuration are not proven optimal, so they get no ground-state claim.

One CSV row per trial goes to standard output; a summary goes to standard error.
"""

import argparse
import csv
import json
from pathlib import Path
import sys

PROVEN_CATEGORIES = {"feasible_fold_minimum", "complete_qubo_global_minimum"}
TOLERANCE = 1e-6  # the same slack the other Python checks use when comparing energies
COLUMNS = (
    "input", "benchmark", "trial", "seed", "reference_category", "reference_status",
    "ground_state_energy", "valid_found", "first_valid_energy", "first_valid_sweep",
    "first_valid_attempted_updates", "first_valid_seconds", "best_valid_energy",
    "ground_state_found", "ground_state_sweep", "ground_state_attempted_updates",
    "ground_state_seconds", "below_reference_investigate",
)


def read_reference(path):
    """Return (reference metadata, variable count) from the first four lines of a QUBO1 file."""
    with path.open(encoding="utf-8") as stream:
        header = [stream.readline().rstrip("\r\n") for _ in range(4)]
    if header[0] != "QUBO1" or not header[3].startswith("meta "):
        raise ValueError(f"{path}: not a QUBO1 file")
    return json.loads(header[3][5:]).get("reference", {}), int(header[1].split()[1])


def blank(value):
    """CSV cell for a value that may not exist."""
    return "" if value is None else value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("logs", nargs="+", type=Path, help="JSONL logs written with --fold-metrics")
    args = parser.parse_args()
    references = {}
    writer = csv.writer(sys.stdout)
    writer.writerow(COLUMNS)
    trials = valid_trials = proven_trials = ground_trials = 0
    for path in args.logs:
        with path.open(encoding="utf-8") as stream:
            for line_number, line in enumerate(stream, 1):
                if not line.strip():
                    continue
                record = json.loads(line)
                if "valid_history" not in record:
                    parser.error(f"{path}:{line_number}: no fold metrics; rerun with --fold-metrics")
                source = Path(record["input"])
                if source not in references:
                    try:
                        references[source] = read_reference(source)
                    except (OSError, ValueError) as error:
                        parser.error(f"{path}:{line_number}: cannot read the input file: {error}")
                reference, variables = references[source]
                if variables != record["variables"]:
                    parser.error(f"{path}:{line_number}: {source} no longer matches this log")

                history = record["valid_history"]  # oldest first, energies strictly decreasing
                first = history[0] if history else None
                proven = reference.get("category") in PROVEN_CATEGORIES
                ground = reference["energy"] if proven else None
                found = first_hit = below = None
                if proven:
                    reached = [entry for entry in history if entry["energy"] <= ground + TOLERANCE]
                    found = bool(reached)
                    first_hit = reached[0] if reached else None  # first time the minimum was reached
                    below = any(entry["energy"] < ground - TOLERANCE for entry in history)
                writer.writerow((
                    record["input"], record["benchmark"], record["trial"], record["seed"],
                    blank(reference.get("category")), blank(reference.get("status")), blank(ground),
                    bool(history),
                    blank(first and first["energy"]), blank(first and first["sweep"]),
                    blank(first and first["attempted_updates"]), blank(first and first["elapsed_seconds"]),
                    blank(history[-1]["energy"] if history else None),
                    blank(found), blank(first_hit and first_hit["sweep"]),
                    blank(first_hit and first_hit["attempted_updates"]),
                    blank(first_hit and first_hit["elapsed_seconds"]), blank(below),
                ))
                trials += 1
                valid_trials += bool(history)
                proven_trials += proven
                ground_trials += bool(found)
    print(f"{trials} trials; a valid fold was found in {valid_trials}; the proven ground state "
          f"was reached in {ground_trials} of {proven_trials} trials with a proven reference",
          file=sys.stderr)


if __name__ == "__main__":
    main()
