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
#
# Example:
#   bash run_batch.sh 8
#   PIPELINE_TIMEOUT=1800 bash run_batch.sh 8   # tighter 30-min cap

N_JOBS=${1:-8}
PIPELINE_TIMEOUT=${PIPELINE_TIMEOUT:-3600}
BENCHMARK_TSV="$HOME/motif/benchmark/mcsa_representatives_parsed_monomers.tsv"

if [ ! -f "$BENCHMARK_TSV" ]; then
    echo "ERROR: Benchmark list not found at $BENCHMARK_TSV"
    echo "Generate it first with extract_pdb_list_from_mcsa.py"
    exit 1
fi

TOTAL=$(tail -n +2 "$BENCHMARK_TSV" | wc -l)

echo "============================================"
echo "Full Pipeline Run"
echo "  Input:    $BENCHMARK_TSV ($TOTAL entries)"
echo "  Jobs:     $N_JOBS parallel"
echo "  Timeout:  ${PIPELINE_TIMEOUT}s per PDB"
echo "  Est time: ~$((TOTAL * 108 / N_JOBS / 60)) minutes"
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
