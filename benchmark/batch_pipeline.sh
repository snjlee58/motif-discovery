#!/bin/bash
set -euo pipefail

# Batch benchmarking pipeline
# Runs family_pipeline.sh for each M-CSA reference entry and collects results.
#
# Usage:
#   bash batch_pipeline.sh <benchmark_list.tsv> [--resume]
#
# Generate the input list with:
#   python3 generate_benchmark_list.py catalytic_residues_homologues_parsed.tsv -o benchmark_list.tsv
#
# For a pilot run:
#   python3 generate_benchmark_list.py catalytic_residues_homologues_parsed.tsv -o pilot_20.tsv --max 20 --shuffle

TSV=${1:?"Usage: bash batch_pipeline.sh <benchmark_list.tsv>"}
RESUME=${2:-""}  # Pass --resume to skip already-completed entries

MOTIF_DIR="/home/sunnylee/motif"
RESULTS_DIR="$SCRATCH/batch_results"
SUMMARY_FILE="$RESULTS_DIR/batch_summary.tsv"

mkdir -p "$RESULTS_DIR"

# Count total entries
TOTAL=$(tail -n +2 "$TSV" | wc -l)
CURRENT=0
SUCCEEDED=0
FAILED=0
SKIPPED=0

echo "============================================"
echo "Batch Benchmark Pipeline"
echo "  Input:   $TSV ($TOTAL entries)"
echo "  Results: $RESULTS_DIR/"
if [ "$RESUME" = "--resume" ]; then
    echo "  Mode:    RESUME (skipping completed)"
fi
echo "============================================"

# Initialize summary file (added elapsed_sec column)
if [ ! -f "$SUMMARY_FILE" ] || [ "$RESUME" != "--resume" ]; then
    echo -e "mcsa_id\tpdb_id\tn_catalytic_residues\tstatus\tprecision\trecall\tf1\tn_predicted\tn_true\telapsed_sec\terror" > "$SUMMARY_FILE"
fi

tail -n +2 "$TSV" | while IFS=$'\t' read -r MCSA_ID PDB_ID N_RES; do
    CURRENT=$((CURRENT + 1))
    PDB_UPPER=$(echo "$PDB_ID" | tr '[:lower:]' '[:upper:]')

    # Output dir and log file live together
    OUTDIR="batch_family_${PDB_UPPER}"
    mkdir -p "$SCRATCH/$OUTDIR"
    LOG_FILE="$SCRATCH/$OUTDIR/pipeline.log"
    PERF_FILE="$SCRATCH/$OUTDIR/baseline_performance.json"

    echo ""
    echo "[$CURRENT/$TOTAL] $PDB_UPPER (mcsa_id=$MCSA_ID, n_cat_res=$N_RES)"

    # Resume mode: skip if performance file already exists
    if [ "$RESUME" = "--resume" ] && [ -f "$PERF_FILE" ]; then
        echo "  ↩ Already completed, skipping"
        SKIPPED=$((SKIPPED + 1))

        # Still extract metrics for summary
        python3 -c "
import json, sys
with open('$PERF_FILE') as f:
    d = json.load(f)
m = d['metrics']
print(f\"$MCSA_ID\t$PDB_UPPER\t$N_RES\tSUCCESS\t{m['precision']:.4f}\t{m['recall']:.4f}\t{m['f1']:.4f}\t{m['n_predicted']}\t{m['n_true']}\t\t\")
" >> "$SUMMARY_FILE" 2>/dev/null || true
        continue
    fi

    # Download PDB if not already present
    PDB_FILE="$SCRATCH/${PDB_UPPER}.pdb"
    if [ ! -f "$PDB_FILE" ]; then
        echo "  Downloading $PDB_UPPER from RCSB..."
        if ! wget -q --timeout=30 \
            "https://files.rcsb.org/download/${PDB_UPPER}.pdb" \
            -O "$PDB_FILE" 2>/dev/null; then
            echo "  ✗ Could not download PDB, skipping"
            echo -e "$MCSA_ID\t$PDB_UPPER\t$N_RES\tFAILED\t\t\t\t\t\t\tPDB download failed" >> "$SUMMARY_FILE"
            FAILED=$((FAILED + 1))
            rm -f "$PDB_FILE"
            continue
        fi
    fi

    # Run pipeline with timing
    echo "  Running family_pipeline.sh..."
    START_TIME=$(date +%s)

    if (cd "$MOTIF_DIR" && bash family_pipeline.sh "$PDB_UPPER" "" "$OUTDIR" --quiet); then
        END_TIME=$(date +%s)
        ELAPSED=$((END_TIME - START_TIME))
        echo "  ✓ Pipeline completed (${ELAPSED}s)"

        # Extract metrics from baseline_performance.json
        if [ -f "$PERF_FILE" ]; then
            python3 -c "
import json, sys
with open('$PERF_FILE') as f:
    d = json.load(f)
m = d['metrics']
print(f\"$MCSA_ID\t$PDB_UPPER\t$N_RES\tSUCCESS\t{m['precision']:.4f}\t{m['recall']:.4f}\t{m['f1']:.4f}\t{m['n_predicted']}\t{m['n_true']}\t$ELAPSED\t\")
" >> "$SUMMARY_FILE"
            SUCCEEDED=$((SUCCEEDED + 1))
        else
            echo "  ⚠ Pipeline ran but no performance file found"
            echo -e "$MCSA_ID\t$PDB_UPPER\t$N_RES\tNO_RESULT\t\t\t\t\t\t${ELAPSED}\tNo baseline_performance.json" >> "$SUMMARY_FILE"
            FAILED=$((FAILED + 1))
        fi
    else
        END_TIME=$(date +%s)
        ELAPSED=$((END_TIME - START_TIME))
        echo "  ✗ Pipeline failed in ${ELAPSED}s (check $LOG_FILE)"
        # Try to capture the error
        ERROR_MSG=$(tail -1 "$LOG_FILE" 2>/dev/null | tr '\t' ' ' | head -c 100)
        echo -e "$MCSA_ID\t$PDB_UPPER\t$N_RES\tFAILED\t\t\t\t\t\t${ELAPSED}\t${ERROR_MSG}" >> "$SUMMARY_FILE"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "============================================"
echo "Batch Complete"
echo "  Succeeded: $SUCCEEDED"
echo "  Failed:    $FAILED"
echo "  Skipped:   $SKIPPED"
echo "  Summary:   $SUMMARY_FILE"
echo "============================================"
echo ""
echo "Next: python3 aggregate_benchmark.py $SUMMARY_FILE"