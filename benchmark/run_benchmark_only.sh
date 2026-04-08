#!/bin/bash
set -euo pipefail

# Re-benchmark scoring only — no downloads, no FoldMason, no conservation.
# Runs benchmark_mcsa.py on all existing batch_family_* folders in parallel.
# Creates a timestamped results folder with per-entry results and a summary.
#
# Usage:
#   bash run_benchmark_only.sh [n_jobs]
#
# Example:
#   bash run_benchmark_only.sh 8

N_JOBS=${1:-8}
MOTIF_DIR="$HOME/motif"
MCSA_FILE="$SCRATCH/m-csa/catalytic_residues_homologues_parsed.tsv"
TIMESTAMP=$(date +%y%m%d_%H%M%S)
RESULTS_DIR="$SCRATCH/benchmark_results_${TIMESTAMP}"

mkdir -p "$RESULTS_DIR"

if [ ! -f "$MCSA_FILE" ]; then
    echo "ERROR: M-CSA file not found at $MCSA_FILE"
    exit 1
fi

# Find all batch_family dirs that have conservation data
DIRS=()
for dir in $SCRATCH/batch_family_*/; do
    PDB_ID=$(basename "$dir" | sed 's/batch_family_//')
    PDB_LOWER=$(echo "$PDB_ID" | tr '[:upper:]' '[:lower:]')
    CONS="$dir/${PDB_LOWER}_conservation.json"
    MAP="$dir/alignment_mapping.json"
    if [ -f "$CONS" ] && [ -f "$MAP" ]; then
        DIRS+=("$PDB_ID")
    fi
done

TOTAL=${#DIRS[@]}

echo "============================================"
echo "Benchmark Scoring (re-score only)"
echo "  Entries:  $TOTAL"
echo "  Jobs:     $N_JOBS parallel"
echo "  Output:   $RESULTS_DIR/"
echo "  M-CSA:    $MCSA_FILE"
echo "============================================"

# Write PDB IDs to a temp file for parallel
TMPFILE=$(mktemp)
printf '%s\n' "${DIRS[@]}" > "$TMPFILE"

# Run benchmark_mcsa.py in parallel
cat "$TMPFILE" | parallel --progress -j "$N_JOBS" '
    PDB_ID={}
    PDB_LOWER=$(echo "$PDB_ID" | tr "[:upper:]" "[:lower:]")
    DIR="'"$SCRATCH"'/batch_family_${PDB_ID}"
    
    CONS="$DIR/${PDB_LOWER}_conservation.json"
    MAP="$DIR/alignment_mapping.json"
    P2RANK="$DIR/p2rank_scores.json"
    PDB_FILE="'"$SCRATCH"'/pdb_files/${PDB_ID}.pdb"
    MCSA="'"$MCSA_FILE"'"
    OUTDIR="'"$RESULTS_DIR"'"
    
    P2RANK_ARG=""
    [ -f "$P2RANK" ] && P2RANK_ARG="--p2rank-json $P2RANK"
    
    PDB_ARG=""
    [ -f "$PDB_FILE" ] && PDB_ARG="--pdb-file $PDB_FILE"
    
    cd '"$MOTIF_DIR"' && python3 benchmark_mcsa.py \
        "$CONS" "$MCSA" "$MAP" \
        --pdb-id "$PDB_LOWER" \
        --top-n auto \
        --exclude-gaps \
        --catalytic-propensity \
        $P2RANK_ARG \
        $PDB_ARG \
        --min-identity 0.2 \
        --output "$OUTDIR/${PDB_ID}.json" \
        2>/dev/null
'

rm -f "$TMPFILE"

# Collect summary
SUMMARY="$RESULTS_DIR/benchmark_summary.tsv"
echo -e "pdb_id\tprecision\trecall\tf1\tn_predicted\tn_true" > "$SUMMARY"

SUCCEEDED=0
FAILED=0
TOTAL_F1=0

for json_file in "$RESULTS_DIR"/*.json; do
    [ -f "$json_file" ] || continue
    [ "$(basename "$json_file")" = "benchmark_summary.json" ] && continue
    
    PDB_ID=$(basename "$json_file" .json)
    
    RESULT=$(python3 -c "
import json
with open('$json_file') as f:
    d = json.load(f)
m = d['metrics']
print(f\"{m['precision']:.4f}\t{m['recall']:.4f}\t{m['f1']:.4f}\t{m['n_predicted']}\t{m['n_true']}\")
" 2>/dev/null) || true
    
    if [ -n "$RESULT" ]; then
        echo -e "$PDB_ID\t$RESULT" >> "$SUMMARY"
        F1=$(echo "$RESULT" | cut -f3)
        TOTAL_F1=$(python3 -c "print($TOTAL_F1 + $F1)")
        SUCCEEDED=$((SUCCEEDED + 1))
    else
        FAILED=$((FAILED + 1))
    fi
done

# Compute aggregate stats
python3 -c "
import csv, statistics

f1_scores = []
zeros = 0
with open('$SUMMARY') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        f1 = float(row['f1'])
        f1_scores.append(f1)
        if f1 == 0:
            zeros += 1

if f1_scores:
    print()
    print('============================================')
    print('BENCHMARK SUMMARY')
    print('============================================')
    print(f'  Entries scored:  {len(f1_scores)}')
    print(f'  Mean F1:         {statistics.mean(f1_scores):.4f}')
    print(f'  Median F1:       {statistics.median(f1_scores):.4f}')
    print(f'  Std F1:          {statistics.stdev(f1_scores):.4f}' if len(f1_scores) > 1 else '')
    print(f'  F1 = 0:          {zeros} ({100*zeros/len(f1_scores):.1f}%)')
    print(f'  F1 >= 0.5:       {sum(1 for x in f1_scores if x >= 0.5)} ({100*sum(1 for x in f1_scores if x >= 0.5)/len(f1_scores):.1f}%)')
    print(f'  F1 >= 0.8:       {sum(1 for x in f1_scores if x >= 0.8)} ({100*sum(1 for x in f1_scores if x >= 0.8)/len(f1_scores):.1f}%)')
    print(f'  Perfect (1.0):   {sum(1 for x in f1_scores if x == 1.0)}')
    print(f'============================================')
    print(f'  Results:  $RESULTS_DIR/')
    print(f'  Summary:  $SUMMARY')
    print(f'============================================')
"

echo ""
echo "Failed: $FAILED"
