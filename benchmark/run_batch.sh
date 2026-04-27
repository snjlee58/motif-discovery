#!/bin/bash
set -euo pipefail

# Run full pipeline on all representative monomers in parallel.
# This runs pipeline.sh for each entry (downloads, FoldMason, conservation, benchmark).
# Only need to rerun this when the pipeline itself changes (not for scoring tweaks).
#
# Usage:
#   bash run_batch.sh [n_jobs]
#
# Tunables (env vars):
#   PIPELINE_TIMEOUT  per-pipeline wallclock cap in seconds (default: 3600 = 1h).
#                     A single stuck PDB (huge cluster, hung download) gets killed
#                     after this so it can't block its xargs slot indefinitely.
#   SKIP_EXISTING     path to an existing results dir. PDBs that already have
#                     <SKIP_EXISTING>/<PDB>/baseline_performance.json are filtered
#                     out so we only re-run the missing/failed ones.
#
# Examples:
#   bash run_batch.sh 8
#   PIPELINE_TIMEOUT=1800 bash run_batch.sh 8   # tighter 30-min cap
#   SKIP_EXISTING=/fast/sunny/motif/batch_runs/<old_dir> bash run_batch.sh 8

N_JOBS=${1:-8}
PIPELINE_TIMEOUT=${PIPELINE_TIMEOUT:-3600}
SKIP_EXISTING="${SKIP_EXISTING:-}"
BENCHMARK_TSV="$HOME/motif/benchmark/mcsa_representatives_parsed_monomers.tsv"

if [ ! -f "$BENCHMARK_TSV" ]; then
    echo "ERROR: Benchmark list not found at $BENCHMARK_TSV"
    echo "Generate it first with extract_pdb_list_from_mcsa.py"
    exit 1
fi

# If SKIP_EXISTING is set, drop PDBs that already have a baseline_performance.json
# in that dir. Lets us retry just the failed/missing PDBs from a prior run.
if [ -n "$SKIP_EXISTING" ]; then
    [ ! -d "$SKIP_EXISTING" ] && { echo "ERROR: SKIP_EXISTING dir not found: $SKIP_EXISTING"; exit 1; }
    FILTERED_TSV=$(mktemp)
    head -1 "$BENCHMARK_TSV" > "$FILTERED_TSV"
    while IFS=$'\t' read -r mcsa pdb nres; do
        [ ! -f "$SKIP_EXISTING/$pdb/baseline_performance.json" ] && \
            printf '%s\t%s\t%s\n' "$mcsa" "$pdb" "$nres" >> "$FILTERED_TSV"
    done < <(tail -n +2 "$BENCHMARK_TSV")
    BEFORE=$(tail -n +2 "$BENCHMARK_TSV" | wc -l)
    AFTER=$(tail -n +2 "$FILTERED_TSV" | wc -l)
    echo "SKIP_EXISTING: $((BEFORE - AFTER)) already done in $SKIP_EXISTING, retrying $AFTER"
    BENCHMARK_TSV="$FILTERED_TSV"
fi

TOTAL=$(tail -n +2 "$BENCHMARK_TSV" | wc -l)

echo "============================================"
echo "Full Pipeline Run"
echo "  Input:    $BENCHMARK_TSV ($TOTAL entries)"
echo "  Jobs:     $N_JOBS parallel"
echo "  Timeout:  ${PIPELINE_TIMEOUT}s per PDB"
echo "  Est time: ~$((TOTAL * 108 / N_JOBS / 60)) minutes"
[ -n "$SKIP_EXISTING" ] && echo "  Retry mode: skipping PDBs already in $SKIP_EXISTING"
echo "============================================"

# Extract column 2 (PDB ID) and run pipeline.sh on each, N_JOBS at a time.
# Using xargs -P rather than GNU parallel so it works without extra deps.
# Per-PDB outputs land in $SCRATCH/batch/<PDB>/ — the `batch/` parent gives
# the persist trap a clean glob target.
# `timeout` sends SIGTERM after PIPELINE_TIMEOUT seconds; pipeline.sh dies,
# xargs moves on to the next PDB. Partial results from killed runs may still
# get persisted by the sbatch trap (whatever made it to disk).
export PIPELINE_TIMEOUT
tail -n +2 "$BENCHMARK_TSV" | cut -f2 | \
    xargs -P "$N_JOBS" -I{} bash -c \
    'cd ~/motif && timeout "${PIPELINE_TIMEOUT}s" bash pipeline.sh "$1" "batch/$1" --quiet' _ {}

echo ""
echo "Done! Results in \$SCRATCH/batch/*/"
