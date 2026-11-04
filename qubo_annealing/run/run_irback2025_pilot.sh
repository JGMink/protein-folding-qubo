#!/usr/bin/env bash
# Bounded pilot only: 8 trials, 2,534,400,000 attempted single-bit updates.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
mkdir -p build/audit experiments/irback2025_audit
make -C qubo_annealing
for sequence in mj5 mj1; do
    python3 -m qubo_bench.export --preset "irback2025_${sequence}" \
        --penalties 1.5 2 2 --output "build/audit/${sequence}.qubo"
    for sweeps in 1000 10000; do
        qubo_annealing/qubo_anneal "build/audit/${sequence}.qubo" \
            --schedule paper25 --selection random --sweeps-per-temperature "$sweeps" \
            --trials 2 --seed 42 \
            --output "experiments/irback2025_audit/${sequence}_${sweeps}.jsonl"
    done
done
python3 experiments/irback2025_audit/verify.py
