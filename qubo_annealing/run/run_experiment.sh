#!/usr/bin/env bash
# One experiment; all annealing settings are ordinary solver arguments.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ $# -lt 2 ]]; then
    echo "Usage: bash run_experiment.sh BENCHMARK_OR_QUBO RUN_ID [solver options]" >&2
    exit 2
fi
INPUT="$1"
RUN_ID="$2"
shift 2
case "$RUN_ID" in
    ""|.|..|*/*|*\\*) echo "RUN_ID must be a single directory name" >&2; exit 2 ;;
esac
for argument in "$@"; do
    if [[ "$argument" == --output ]]; then
        echo "Output is managed by the run directory; do not pass --output" >&2
        exit 2
    fi
done
case "$INPUT" in
    fold8|fold24|fold32) INPUT="$ROOT/benchmarks/$INPUT.qubo" ;;
    mj1|mj2|mj3|mj4|mj5|mj6)
        INPUT="${MJ_DIR:-$ROOT/benchmarks/mj48}/$INPUT.qubo" ;;
esac
if [[ ! -f "$INPUT" ]]; then
    echo "Missing QUBO: $INPUT (prepare MJ inputs with prepare_mj48.sh)" >&2
    exit 2
fi
make -C "$ROOT"
RESULTS_DIR="${RESULTS_DIR:-$ROOT/results}"
mkdir -p "$RESULTS_DIR"
RUN_DIR="$RESULTS_DIR/$RUN_ID"
# mkdir fails if the identifier already exists, including an incomplete prior run.
mkdir "$RUN_DIR"
sha256sum "$INPUT" "$ROOT/qubo_anneal" > "$RUN_DIR/checksums.sha256"
cp "$ROOT/Makefile" "$RUN_DIR/Makefile.snapshot"
printf '%s
' "$INPUT" "$@" > "$RUN_DIR/arguments.txt"
date -u > "$RUN_DIR/started_utc.txt"
"$ROOT/qubo_anneal" "$INPUT" "$@" --output "$RUN_DIR/trials.jsonl"
printf 'complete
' > "$RUN_DIR/status.txt"
echo "Results: $RUN_DIR/trials.jsonl"
