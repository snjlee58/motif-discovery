#!/bin/bash
#SBATCH --job-name=motif-single
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=/home/sunnylee/motif/slurm_logs/motif_single_%j.out

# Single-PDB motif-discovery run.
#
# $SCRATCH = /mnt/scratch/$USER/$SLURM_JOB_ID is node-local and AUTO-DELETED on
# job exit. The trap copies results to /fast/sunny/motif/single_runs/ before
# Slurm wipes scratch.
#
# Usage: sbatch scripts/single_run.sbatch <PDB_ID>
# Example: sbatch scripts/single_run.sbatch 1N8I

set -uo pipefail

PDB_ID=${1:?"Usage: sbatch scripts/single_run.sbatch <PDB_ID>"}

# sbatch runs a bare non-interactive shell — load conda's shell functions
# explicitly so `conda activate` is defined.
source ~/miniforge3/etc/profile.d/conda.sh
conda activate motif

export FAST="/fast/sunny"

TS=$(date +%y%m%d_%H%M%S)
PERSIST="/fast/sunny/motif/single_runs/${TS}_job${SLURM_JOB_ID}_${PDB_ID}"
mkdir -p "$PERSIST"

persist_results() {
    echo ""
    echo "=================================="
    echo "Persisting results to $PERSIST..."
    cp -r "$SCRATCH"/*_family_"$PDB_ID" "$PERSIST"/ 2>/dev/null || true
    n=$(ls -d "$PERSIST"/*_family_"$PDB_ID"/ 2>/dev/null | wc -l)
    echo "Persisted $n result dir(s) to $PERSIST"
    echo "Persistence done: $(date)"
    echo "=================================="
}
trap persist_results EXIT

echo "=================================="
echo "Job ID:     $SLURM_JOB_ID"
echo "Node:       $(hostname)"
echo "PDB ID:     $PDB_ID"
echo "CPUs:       $SLURM_CPUS_PER_TASK"
echo "SCRATCH:    $SCRATCH  (auto-deleted at job end)"
echo "FAST:       $FAST"
echo "Persist:    $PERSIST"
echo "Started:    $(date)"
echo "=================================="

cd ~/motif
bash pipeline.sh "$PDB_ID"
pipeline_exit=$?

echo "=================================="
echo "pipeline.sh exited with code: $pipeline_exit"
echo "Pipeline finished: $(date)"
