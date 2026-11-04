#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON="${PYTHON:-python3}"
make -C "$ROOT"
mkdir -p "$ROOT/results"
cd "$ROOT"
./qubo_anneal benchmarks/fold8.qubo --sweeps 2 --trials 1 --seed 7 \
    --checkpoint 1 --target -2.84 --output results/smoke.jsonl
"$PYTHON" - <<'PY'
import json
from pathlib import Path

records = [json.loads(line) for line in Path("results/smoke.jsonl").read_text().splitlines()]
assert len(records) == 1
record = records[0]
assert record["benchmark"] == "fold8" and record["attempted_updates"] == 128
assert len(record["best_bits"]) == len(record["final_bits"]) == 64
print("Smoke run passed: results/smoke.jsonl")
PY
