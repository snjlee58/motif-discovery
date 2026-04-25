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
# $SCRATCH under SLURM = /mnt/scratch/$USER/$SLURM_JOB_ID and is AUTO-DELETED
# on job exit. We persist small result files to $FAST before exiting so they
# survive. The trap fires on normal exit, Ctrl+C, or SIGTERM (within SLURM's
# grace period before SIGKILL).
#
# Submit with: sbatch scripts/motif_benchmark.sbatch

# Note: no -e — we want persistence to run even if the pipeline fails.
set -uo pipefail

# Activate the conda env that has biopython/numpy/etc.
# sbatch runs in a bare non-interactive shell, so conda's shell functions
# aren't auto-sourced from .bashrc. Source them explicitly.
source ~/miniforge3/etc/profile.d/conda.sh
conda activate motif

export FAST="/fast/sunny"
TS=$(date +%y%m%d_%H%M%S)
BATCH_OUT="/fast/sunny/motif/batch_runs/${TS}_job${SLURM_JOB_ID}"
mkdir -p "$BATCH_OUT"

persist_results() {
    echo ""
    echo "=================================="
    echo "Persisting results to $BATCH_OUT..."
    # Portable cp-based copy — no rsync needed. Grabs the small result files
    # (json / fa / nw / txt / log) and skips large intermediate dirs
    # (structures/, p2rank_output/, tmp/) which are re-derivable.
    for d in "$SCRATCH"/batch/*/; do
        [ -d "$d" ] || continue
        name=$(basename "$d")
        dst="$BATCH_OUT/$name"
        mkdir -p "$dst"
        cp "$d"/*.json "$d"/*.fa "$d"/*.nw "$d"/*.txt "$d"/*.log "$dst/" 2>/dev/null || true
    done
    n=$(ls -d "$BATCH_OUT"/*/ 2>/dev/null | wc -l)
    echo "Persisted $n directories to $BATCH_OUT"
    echo "Persistence done: $(date)"
    echo "=================================="
    echo "Aggregate with:"
    echo "  python3 ~/motif/benchmark/summarize_results.py \\"
    echo "    ~/motif/benchmark/mcsa_representatives_parsed_monomers.tsv \\"
    echo "    --scratch $BATCH_OUT"
    echo "=================================="
}
trap persist_results EXIT

echo "=================================="
echo "Job ID:     $SLURM_JOB_ID"
echo "Node:       $(hostname)"
echo "CPUs:       $SLURM_CPUS_PER_TASK"
echo "SCRATCH:    $SCRATCH  (auto-deleted at job end)"
echo "FAST:       $FAST"
echo "Persist:    $BATCH_OUT"
echo "Started:    $(date)"
echo "=================================="

echo "Dependency check on $(hostname):"
for cmd in wget foldmason python3 xargs cp mkdir; do
    if command -v "$cmd" &>/dev/null; then
        echo "  OK   $cmd: $(command -v $cmd)"
    else
        echo "  MISS $cmd: NOT FOUND"
    fi
done

echo ""
echo "Python module check:"
for mod in Bio numpy; do
    if python3 -c "import $mod" &>/dev/null; then
        echo "  OK   python3 -c 'import $mod'"
    else
        echo "  MISS python3 cannot import $mod  (run: pip install --user biopython numpy)"
    fi
done

echo ""
echo "Network check (AFDB cluster file):"
if [ -f "$FAST/afdb_clusters/5-allmembers-repId-entryId-cluFlag-taxId.tsv" ]; then
    echo "  OK   cluster file present"
else
    echo "  MISS cluster file not found at $FAST/afdb_clusters/"
fi

echo ""
echo "Network check (UniProt API reachable):"
if wget -q --timeout=10 --spider https://rest.uniprot.org/idmapping/run 2>/dev/null; then
    echo "  OK   rest.uniprot.org reachable"
else
    # spider might fail with POST-only endpoints; try a simpler URL
    if wget -q --timeout=10 -O /dev/null https://rest.uniprot.org/ 2>/dev/null; then
        echo "  OK   rest.uniprot.org reachable"
    else
        echo "  WARN rest.uniprot.org NOT reachable from this node"
    fi
fi

echo ""
echo "Network check (AlphaFold EBI reachable):"
if wget -q --timeout=10 -O /dev/null https://alphafold.ebi.ac.uk/ 2>/dev/null; then
    echo "  OK   alphafold.ebi.ac.uk reachable"
else
    echo "  WARN alphafold.ebi.ac.uk NOT reachable from this node"
fi
echo "=================================="

cd ~/motif/benchmark
bash run_batch.sh "$SLURM_CPUS_PER_TASK"
pipeline_exit=$?

echo "=================================="
echo "run_batch.sh exited with code: $pipeline_exit"
echo "Pipeline finished: $(date)"
