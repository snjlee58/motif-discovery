#!/bin/bash
#SBATCH --job-name=motif-benchmark
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/home/sunnylee/motif/slurm_logs/motif_%j.out

# Run the full motif-discovery benchmark across all M-CSA representative monomers.
#
# IMPORTANT: per the lab guide, $SCRATCH under SLURM resolves to
# /mnt/scratch/$USER/$SLURM_JOB_ID and is AUTO-DELETED when the job ends.
# We use it for heavy I/O (downloads, FoldMason MSAs) but rsync the
# non-reproducible outputs to $FAST before exit so they survive.
#
# Submit with:  sbatch scripts/motif_benchmark.sbatch
# Monitor with: squeue -u $USER  /  tail -f slurm_logs/motif_<jobid>.out

set -euo pipefail

export FAST="/fast/sunny"
OUTPUT_PERSIST="/fast/sunny/batch_out/job_${SLURM_JOB_ID}"
mkdir -p "$OUTPUT_PERSIST"

# Rsync results to persistent storage before $SCRATCH is auto-deleted.
# Runs on normal exit, Ctrl+C, or SIGTERM from SLURM's grace period.
persist_results() {
    echo ""
    echo "=================================="
    echo "Persisting results to $OUTPUT_PERSIST..."
    # Skip structures/ (re-downloadable from AFDB) and intermediate p2rank/tmp dirs
    rsync -a \
        --exclude='structures/' \
        --exclude='p2rank_output/' \
        --exclude='tmp/' \
        "$SCRATCH/" "$OUTPUT_PERSIST/" 2>&1 || true
    echo "Persistence done: $(date)"
    echo "=================================="
    echo "To aggregate results, run from anywhere:"
    echo "  SCRATCH=$OUTPUT_PERSIST python3 ~/motif/benchmark/summarize_benchmark.py \\"
    echo "    ~/motif/benchmark/benchmark_representative_monomers.tsv"
    echo "=================================="
}
trap persist_results EXIT

echo "=================================="
echo "Job ID:     $SLURM_JOB_ID"
echo "Node:       $(hostname)"
echo "CPUs:       $SLURM_CPUS_PER_TASK"
echo "SCRATCH:    $SCRATCH  (auto-deleted at job end)"
echo "FAST:       $FAST"
echo "Persist:    $OUTPUT_PERSIST"
echo "Started:    $(date)"
echo "=================================="

cd ~/motif/benchmark
bash run_full_pipeline.sh "$SLURM_CPUS_PER_TASK"

echo "=================================="
echo "Pipeline finished: $(date)"
