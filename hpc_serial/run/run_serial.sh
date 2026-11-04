#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Edit experiment settings here. Each output line also stores these settings.
BENCHMARKS=(fold8 fold24 fold32)
SCHEDULES=(geometric)
SWEEPS_FOLD8=(1000 5000)
SWEEPS_FOLD24=(5000 20000)
SWEEPS_FOLD32=(5000 20000)
TRIALS=5
SEED=42
# fold8 retains the repository's tested range. The larger ranges are rounded
# inverses of D-Wave's automatic beta bounds for these exact fixed QUBOs.
TEMP_START_FOLD8=10
TEMP_END_FOLD8=0.1
TEMP_START_FOLD24=570
TEMP_END_FOLD24=0.00885
TEMP_START_FOLD32=389
TEMP_END_FOLD32=0.00818
CHECKPOINT=0
OUTPUT_DIR="$ROOT/results"
# Empty values disable the optional numeric QUBO energy threshold.
# Reaching a threshold does not establish that the saved fold is feasible.
TARGET_FOLD8=-2.84
TARGET_FOLD24=-10.05
TARGET_FOLD32=-14.25

make -C "$ROOT"
mkdir -p "$OUTPUT_DIR"
cd "$ROOT"

for benchmark in "${BENCHMARKS[@]}"; do
    case "$benchmark" in
        fold8)
            target="$TARGET_FOLD8"; sweeps_list=("${SWEEPS_FOLD8[@]}")
            temp_start="$TEMP_START_FOLD8"; temp_end="$TEMP_END_FOLD8" ;;
        fold24)
            target="$TARGET_FOLD24"; sweeps_list=("${SWEEPS_FOLD24[@]}")
            temp_start="$TEMP_START_FOLD24"; temp_end="$TEMP_END_FOLD24" ;;
        fold32)
            target="$TARGET_FOLD32"; sweeps_list=("${SWEEPS_FOLD32[@]}")
            temp_start="$TEMP_START_FOLD32"; temp_end="$TEMP_END_FOLD32" ;;
        *) echo "Unknown benchmark: $benchmark" >&2; exit 2 ;;
    esac
    for schedule in "${SCHEDULES[@]}"; do
        for sweeps in "${sweeps_list[@]}"; do
            output="$OUTPUT_DIR/${benchmark}_${schedule}_${sweeps}.jsonl"
            if [[ -n "$target" ]]; then
                ./serial "benchmarks/${benchmark}.qubo" \
                    --schedule "$schedule" --temp-start "$temp_start" --temp-end "$temp_end" \
                    --sweeps "$sweeps" --trials "$TRIALS" --seed "$SEED" \
                    --checkpoint "$CHECKPOINT" --target "$target" --output "$output"
            else
                ./serial "benchmarks/${benchmark}.qubo" \
                    --schedule "$schedule" --temp-start "$temp_start" --temp-end "$temp_end" \
                    --sweeps "$sweeps" --trials "$TRIALS" --seed "$SEED" \
                    --checkpoint "$CHECKPOINT" --output "$output"
            fi
            echo "Wrote $output"
        done
    done
done
