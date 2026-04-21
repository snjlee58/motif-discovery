#!/bin/bash
set -euo pipefail

# Run full pipeline on all representative monomers in parallel.
# This runs pipeline.sh for each entry (downloads, FoldMason, conservation, benchmark).
# Only need to rerun this when the pipeline itself changes (not for scoring tweaks).
#
# Usage:
#   bash run_full_pipeline.sh [n_jobs]
#
# Example:
#   bash run_full_pipeline.sh 8

N_JOBS=${1:-8}
BENCHMARK_TSV="$HOME/motif/benchmark/benchmark_representative_monomers.tsv"

if [ ! -f "$BENCHMARK_TSV" ]; then
    echo "ERROR: Benchmark list not found at $BENCHMARK_TSV"
    echo "Generate it first with generate_benchmark_list.py"
    exit 1
fi

TOTAL=$(tail -n +2 "$BENCHMARK_TSV" | wc -l)

echo "============================================"
echo "Full Pipeline Run"
echo "  Input:    $BENCHMARK_TSV ($TOTAL entries)"
echo "  Jobs:     $N_JOBS parallel"
echo "  Est time: ~$((TOTAL * 108 / N_JOBS / 60)) minutes"
echo "============================================"

tail -n +2 "$BENCHMARK_TSV" | \
    parallel --progress -j "$N_JOBS" --colsep '\t' \
    'cd ~/motif && bash pipeline.sh {2} "" batch_family_{2} --quiet'

echo ""
echo "Done! Results in \$SCRATCH/batch_family_*/"
echo "Now run: bash run_benchmark_only.sh"
