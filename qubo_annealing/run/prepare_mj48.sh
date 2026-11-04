#!/usr/bin/env bash
# Generate selected paper inputs; coefficients are fixed here, never in the runner.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PYTHON="${PYTHON:-python3}"
MJ_DIR="${MJ_DIR:-$ROOT/qubo_annealing/benchmarks/mj48}"
if [[ $# -eq 0 ]]; then
    set -- mj1 mj2 mj3 mj4 mj5 mj6
fi
mkdir -p "$MJ_DIR"
MJ_DIR="$(cd "$MJ_DIR" && pwd)"
cd "$ROOT"
for benchmark in "$@"; do
    case "$benchmark" in
        mj1|mj2|mj3|mj4|mj5|mj6) ;;
        *) echo "Unknown MJ benchmark: $benchmark" >&2; exit 2 ;;
    esac
    output="$MJ_DIR/$benchmark.qubo"
    if [[ -e "$output" ]]; then
        echo "Refusing to replace $output; select a fresh MJ_DIR" >&2
        exit 2
    fi
    "$PYTHON" -m qubo_bench.export --preset "irback2025_$benchmark" \
        --penalties 1.5 2 2 --output "$output"
done
