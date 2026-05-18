#!/bin/bash
set -euo pipefail

# Pipeline: query → foldseek search (AFDB50) → top-N homologs → foldmason MSA →
#           conservation + P2Rank + 3Di + spatial-clustering scoring →
#           M-CSA benchmark
#
# Prerequisites:
#   - foldseek binary on PATH
#   - AFDB50 foldseek database at $FOLDSEEK_DB
#     (default: /fast/databases/foldseek/afdb/afdb50)
#   - M-CSA TSV at $MCSA_FILE (for benchmarking)
#
# Tunables (env vars; all have defaults):
#   FOLDSEEK_DB, FOLDSEEK_MAX_SEQS, FOLDSEEK_PROB_MIN, FOLDSEEK_TOP_N
#   SCORE_T (combined-score threshold for catalytic prediction)
#   DOWNLOAD_JOBS (parallel AlphaFold downloads)
#
# The AFDB-cluster lookup mode (query → AFDB cluster → all members → MSA) is
# preserved as commented-out code in STEP 1 for the future "analyse a given
# AFDB cluster" use case.
#
# Usage:
#   bash pipeline.sh <pdb_id> [output_dir] [--quiet]
#
# The UniProt ID is resolved automatically from the PDB ID via the UniProt
# ID Mapping API. To override (e.g. if the API is down or maps wrong):
#   UNIPROT_ID=P62593 bash pipeline.sh 1BTL
#
# Examples:
#   bash pipeline.sh 1BTL
#   bash pipeline.sh 1BTL 250212_family_1btl
#   bash pipeline.sh 1BTL 250212_family_1btl --quiet

PIPELINE_START=$(date +%s)

#####################
# ARGUMENTS
#####################
PDB_ID=${1:?"Usage: bash pipeline.sh <pdb_id> [output_dir] [--quiet]"}
OUTDIR=${2:-$(date +%y%m%d_%H%M%S)_family_${PDB_ID}}
QUIET=${3:-""}
UNIPROT_ID=${UNIPROT_ID:-""}   # env-var override; otherwise resolved below

# Combined-score threshold for catalytic-residue selection. Predicted set =
# all residues with combined_score >= SCORE_T (no ground-truth-anchoring).
# Override via env: SCORE_T=1.2 bash pipeline.sh 1BTL
SCORE_T="${SCORE_T:-1.0}"

PDB_ID_LOWER=$(echo "$PDB_ID" | tr '[:upper:]' '[:lower:]')

# $FAST points to shared NVMe storage (e.g. /fast/sunny) holding persistent reference
# databases (AFDB clusters, M-CSA, PDB cache) — visible from every compute node.
# Falls back to $SCRATCH for single-node setups where nothing has been moved yet.
FAST="${FAST:-$SCRATCH}"

CLUSTER_FILE="$FAST/afdb_clusters/v6/5-allmembers-repId-entryId-cluFlag-taxId.tsv"
MCSA_FILE="$FAST/m-csa/catalytic_residues_homologues_parsed.tsv"
PDB_CACHE="$FAST/pdb_files"
mkdir -p "$PDB_CACHE"

# Foldseek homolog-search settings (default workflow).
# AFDB50 database lives in the shared /fast/databases tree, not per-user $FAST.
FOLDSEEK_DB="${FOLDSEEK_DB:-/fast/databases/foldseek/afdb/afdb50}"
FOLDSEEK_MAX_SEQS="${FOLDSEEK_MAX_SEQS:-1000}"   # foldseek --max-seqs
FOLDSEEK_PROB_MIN="${FOLDSEEK_PROB_MIN:-0.8}"    # parse_foldseek_hits --prob-min
FOLDSEEK_TOP_N="${FOLDSEEK_TOP_N:-200}"          # parse_foldseek_hits --top-n

# Set up logging
mkdir -p $SCRATCH/$OUTDIR
LOG_FILE="$SCRATCH/$OUTDIR/pipeline.log"

if [ "$QUIET" = "--quiet" ]; then
    exec > "$LOG_FILE" 2>&1
else
    exec > >(tee "$LOG_FILE") 2>&1
fi

#####################
# STEP 0: Resolve UniProt ID if not provided
#####################
if [ -z "$UNIPROT_ID" ]; then
  echo "Resolving UniProt ID for PDB: $PDB_ID via UniProt ID Mapping API..."

  # Step 1: Submit mapping job (PDB -> UniProtKB)
  SUBMIT_RESPONSE=$(wget -qO- \
    --post-data="from=PDB&to=UniProtKB&ids=${PDB_ID_LOWER}" \
    https://rest.uniprot.org/idmapping/run)
  echo "  Submit response: $SUBMIT_RESPONSE"

  JOB_ID=$(echo "$SUBMIT_RESPONSE" | python3 -c "import sys,json; print(json.load(sys.stdin)['jobId'])" 2>/dev/null)

  if [ -z "$JOB_ID" ]; then
    echo "ERROR: Failed to submit UniProt ID mapping job"
    exit 1
  fi
  echo "  Submitted job: $JOB_ID"

  # Step 2: Poll until complete (wget follows redirects by default)
  for i in $(seq 1 15); do
    sleep 1
    RESULT_RAW=$(wget -qO- "https://rest.uniprot.org/idmapping/status/${JOB_ID}")

    # Check if we got results directly (status redirects to results when done)
    UNIPROT_ID=$(echo "$RESULT_RAW" | python3 -c "
import sys, json
data = json.load(sys.stdin)
# Case 1: redirected to results (has 'results' key)
results = data.get('results', [])
if results:
    # Could be full UniProtKB entry or simple mapping
    to = results[0].get('to', {})
    if isinstance(to, dict):
        print(to.get('primaryAccession', ''))
    else:
        print(to)
    sys.exit(0)
# Case 2: still running
status = data.get('jobStatus', '')
if status == 'RUNNING':
    print('__RUNNING__')
else:
    print('')
" 2>/dev/null)

    if [ "$UNIPROT_ID" = "__RUNNING__" ]; then
      echo "  Waiting for job to complete... (${i}s)"
      UNIPROT_ID=""
      continue
    elif [ -n "$UNIPROT_ID" ]; then
      break
    fi

    # Fallback: try results/stream endpoint directly
    UNIPROT_ID=$(wget -qO- "https://rest.uniprot.org/idmapping/results/stream/${JOB_ID}" \
      | python3 -c "
import sys, json
data = json.load(sys.stdin)
results = data.get('results', [])
if results:
    to = results[0].get('to', {})
    if isinstance(to, dict):
        print(to.get('primaryAccession', ''))
    else:
        print(to)
else:
    print('')
" 2>/dev/null)

    if [ -n "$UNIPROT_ID" ]; then
      break
    fi
  done

  if [ -z "$UNIPROT_ID" ]; then
    echo "ERROR: Could not resolve UniProt ID for PDB $PDB_ID"
    echo "  Debug: last raw response from API:"
    echo "  $RESULT_RAW" | head -5
    echo "  Try providing it manually: UNIPROT_ID=<uniprot_id> bash pipeline.sh $PDB_ID"
    exit 1
  fi
  echo "  Resolved: $PDB_ID -> $UNIPROT_ID"
fi

echo "============================================"
echo "Homolog-Search Pipeline (foldseek → foldmason)"
echo "  UniProt:        $UNIPROT_ID"
echo "  PDB:            $PDB_ID"
echo "  Foldseek DB:    $FOLDSEEK_DB"
echo "  Output:         $SCRATCH/$OUTDIR"
echo "============================================"

#####################
# STEP 0b: Download experimental query PDB (needed as foldseek query)
#####################
echo ""
echo "[0b] Ensuring experimental query PDB is available..."

PDB_FILE="$PDB_CACHE/${PDB_ID}.pdb"
if [ ! -f "$PDB_FILE" ]; then
  echo "  Downloading $PDB_ID from RCSB..."
  wget -q --timeout=30 \
    "https://files.rcsb.org/download/${PDB_ID}.pdb" \
    -O "$PDB_FILE" 2>/dev/null || rm -f "$PDB_FILE"
fi

# Fallback: if RCSB has no PDB, use the query's AlphaFold model instead.
if [ ! -f "$PDB_FILE" ]; then
  echo "  RCSB has no $PDB_ID; falling back to AlphaFold model for $UNIPROT_ID"
  PDB_FILE="$PDB_CACHE/AF-${UNIPROT_ID}-F1-model_v6.pdb"
  if [ ! -f "$PDB_FILE" ]; then
    wget -q --timeout=30 \
      "https://alphafold.ebi.ac.uk/files/AF-${UNIPROT_ID}-F1-model_v6.pdb" \
      -O "$PDB_FILE" 2>/dev/null || rm -f "$PDB_FILE"
  fi
fi

if [ ! -f "$PDB_FILE" ]; then
  echo "ERROR: Could not obtain query structure (neither $PDB_ID.pdb from RCSB"
  echo "       nor AF-${UNIPROT_ID}-F1-model_v6.pdb from AlphaFold)."
  exit 1
fi
echo "  Query structure: $PDB_FILE"

#####################
# STEP 1: Foldseek homolog search
#####################
echo ""
echo "[1] Running foldseek easy-search against $FOLDSEEK_DB..."

if [ ! -f "${FOLDSEEK_DB}.dbtype" ]; then
  echo "ERROR: Foldseek database not found at $FOLDSEEK_DB"
  echo "       Expected mmseqs2-style files (e.g. ${FOLDSEEK_DB}.dbtype)."
  exit 1
fi

FOLDSEEK_TSV="$SCRATCH/$OUTDIR/foldseek_hits.tsv"
foldseek easy-search \
  "$PDB_FILE" \
  "$FOLDSEEK_DB" \
  "$FOLDSEEK_TSV" \
  "$SCRATCH/tmp_foldseek" \
  --format-output query,target,fident,alnlen,evalue,prob,lddt,alntmscore \
  --max-seqs "$FOLDSEEK_MAX_SEQS"

N_HITS=$(wc -l < "$FOLDSEEK_TSV")
echo "  Foldseek returned $N_HITS hits"

# Parse → filter → top-N → UniProt accession list
python3 src/parse_foldseek_hits.py \
  "$FOLDSEEK_TSV" \
  --prob-min "$FOLDSEEK_PROB_MIN" \
  --top-n "$FOLDSEEK_TOP_N" \
  -o "$SCRATCH/$OUTDIR/homologs.txt"

# ─────────────────────────────────────────────────────────────────────────
# ARCHIVED — AFDB cluster lookup (kept for the "analyse a given AFDB cluster"
# use case; not used in the default homolog-search workflow).
# ─────────────────────────────────────────────────────────────────────────
# if [ ! -f "$CLUSTER_FILE" ]; then
#   echo "ERROR: Cluster file not found at $CLUSTER_FILE"
#   echo "Download it first (FAST=$FAST):"
#   echo "  mkdir -p $FAST/afdb_clusters/v6 && cd $FAST/afdb_clusters/v6"
#   echo "  wget https://afdb-cluster.steineggerlab.workers.dev/v6/5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz"
#   echo "  gunzip 5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz"
#   exit 1
# fi
#
# # File 5 format: repId \t memId \t cluFlag \t taxId
# #   cluFlag=1: AFDB50 representative (Foldseek-clustered) — kept for our MSA
# #   cluFlag=2: non-rep AFDB50 member (≥50% identical near-duplicate) — dropped
# REP_ID=$(grep "${UNIPROT_ID}" "$CLUSTER_FILE" | head -1 | cut -f1) || true
# if [ -z "$REP_ID" ]; then
#   echo "ERROR: $UNIPROT_ID not found in AFDB clusters"
#   exit 1
# fi
# echo "  Cluster representative: $REP_ID"
# COUNTS=$(awk -F'\t' -v rep="$REP_ID" -v out="$SCRATCH/$OUTDIR/cluster_members.txt" '
#   $1 == rep {
#     full++
#     if ($3 == 1) { print $2 > out; afdb50++ }
#   }
#   END { printf "%d %d", afdb50+0, full+0 }
# ' "$CLUSTER_FILE")
# N_AFDB50=$(echo "$COUNTS" | cut -d' ' -f1)
# N_FULL=$(echo "$COUNTS" | cut -d' ' -f2)
# echo "  Cluster members: $N_AFDB50 AFDB50 reps (of $N_FULL total)"

#####################
# STEP 2: Download AlphaFold structures for foldseek homologs
#####################
echo ""
echo "[2] Downloading AlphaFold structures for homologs..."
mkdir -p $SCRATCH/$OUTDIR/structures

TOTAL=$(wc -l < $SCRATCH/$OUTDIR/homologs.txt)

# Build the to-do list (skip already-cached structures from prior runs).
TODO_FILE="$SCRATCH/$OUTDIR/to_download.txt"
> "$TODO_FILE"
while read -r MEMBER_ID; do
  OUTFILE="$SCRATCH/$OUTDIR/structures/AF-${MEMBER_ID}-F1-model_v6.pdb"
  [ -f "$OUTFILE" ] || echo "$MEMBER_ID" >> "$TODO_FILE"
done < $SCRATCH/$OUTDIR/homologs.txt
NEED=$(wc -l < "$TODO_FILE")
CACHED=$((TOTAL - NEED))
echo "  $TOTAL total, $CACHED already cached, $NEED to download"

# Parallel downloads. AFDB EBI handles tens of concurrent connections fine,
# but if the surrounding batch (run_batch.sh) is also parallelizing
# pipelines, lower this via DOWNLOAD_JOBS env var to avoid hammering the API.
DOWNLOAD_JOBS=${DOWNLOAD_JOBS:-8}

download_af_structure() {
  local MEMBER_ID="$1"
  local OUTFILE="$STRUCT_DIR/AF-${MEMBER_ID}-F1-model_v6.pdb"
  local URL="https://alphafold.ebi.ac.uk/files/AF-${MEMBER_ID}-F1-model_v6.pdb"
  wget -q --timeout=10 -O "$OUTFILE" "$URL" 2>/dev/null || rm -f "$OUTFILE"
}
export -f download_af_structure

if [ "$NEED" -gt 0 ]; then
  echo "  Downloading with $DOWNLOAD_JOBS parallel workers..."
  STRUCT_DIR="$SCRATCH/$OUTDIR/structures"
  export STRUCT_DIR

  xargs -a "$TODO_FILE" -P "$DOWNLOAD_JOBS" -I{} bash -c 'download_af_structure "$@"' _ {}
fi

DOWNLOADED=$(ls "$SCRATCH/$OUTDIR/structures/"AF-*.pdb 2>/dev/null | wc -l)
FAILED=$((TOTAL - DOWNLOADED))
rm -f "$TODO_FILE"

echo "  Done: downloaded=$DOWNLOADED  failed=$FAILED"

# Add the query structure (downloaded earlier in STEP 0b) into the MSA pool.
echo "  Copying query structure into structures/"
cp "$PDB_FILE" $SCRATCH/$OUTDIR/structures/

echo "  DEBUG: counting structures"
N_STRUCTURES=$(find $SCRATCH/$OUTDIR/structures -name "*.pdb" -o -name "*.cif" | wc -l)
echo "  Total structures for alignment: $N_STRUCTURES"

if [ "$N_STRUCTURES" -lt 3 ]; then
  echo "ERROR: Too few structures downloaded ($N_STRUCTURES). Need at least 3."
  exit 1
fi

#####################
# STEP 3: FoldMason structural MSA
#####################
echo ""
echo "[3] Running FoldMason..."

foldmason easy-msa \
  $SCRATCH/$OUTDIR/structures \
  $SCRATCH/$OUTDIR/foldmason_result \
  $SCRATCH/tmp

MSA_FILE="$SCRATCH/$OUTDIR/foldmason_result_aa.fa"
echo "  MSA file: $MSA_FILE"

# Quick stats
N_SEQS=$(grep -c "^>" "$MSA_FILE")
ALN_LEN=$(head -2 "$MSA_FILE" | tail -1 | tr -d '\n' | wc -c)
echo "  Sequences: $N_SEQS"
echo "  Alignment length: $ALN_LEN"

#####################
# STEP 4: Conservation scoring
#####################
echo ""
echo "[4] Running conservation scoring..."
python3 src/score_conservation.py $SCRATCH/$OUTDIR $MSA_FILE $PDB_ID_LOWER

#####################
# STEP 5: Create alignment mapping
#####################
echo ""
echo "[5] Mapping alignment columns to PDB residue IDs..."
python3 src/map_alignment_to_pdb.py \
  $MSA_FILE \
  $PDB_ID \
  --pdb-file $PDB_CACHE/${PDB_ID}.pdb \
  --uniprot $UNIPROT_ID \
  -o $SCRATCH/$OUTDIR/alignment_mapping.json

#####################
# STEP 5b: P2Rank binding site prediction
#####################
echo ""
echo "[5b] Running P2Rank binding site prediction..."

P2RANK_JSON="$SCRATCH/$OUTDIR/p2rank_scores.json"
PDB_FILE="$PDB_CACHE/${PDB_ID}.pdb"

# P2Rank is optional — any failure (missing binary, bad PDB, parser hiccup) must
# not kill the pipeline. Use a function that returns non-zero on any failure,
# then treat that as "P2Rank unavailable" and continue without binding-site data.
run_p2rank() {
  if ! command -v prank &>/dev/null; then
    echo "  SKIP: prank not in PATH"
    return 1
  fi
  if [ ! -f "$PDB_FILE" ]; then
    echo "  SKIP: PDB file missing at $PDB_FILE"
    return 1
  fi
  if ! prank predict -f "$PDB_FILE" -o "$SCRATCH/$OUTDIR/p2rank_output"; then
    echo "  WARNING: prank predict failed"
    return 1
  fi
  local csv
  csv=$(find "$SCRATCH/$OUTDIR/p2rank_output" -name "*_residues.csv" 2>/dev/null | head -1)
  if [ -z "$csv" ]; then
    echo "  WARNING: prank ran but produced no residues CSV"
    return 1
  fi
  if ! python3 src/parse_p2rank.py "$csv" -o "$P2RANK_JSON"; then
    echo "  WARNING: parse_p2rank.py failed"
    return 1
  fi
  echo "  ✓ P2Rank scores saved to: p2rank_scores.json"
  return 0
}

if ! run_p2rank; then
  P2RANK_JSON=""
fi

#####################
# STEP 6: Extract top conserved positions
#####################
echo ""
echo "[6] Extracting top conserved positions..."

# Top 5 with mapping
python3 src/extract_top_conserved.py \
  $SCRATCH/$OUTDIR/${PDB_ID_LOWER}_conservation.json \
  --top-n 5 \
  --exclude-gaps \
  --min-identity 0.2 \
  --mapping $SCRATCH/$OUTDIR/alignment_mapping.json \
  --output $SCRATCH/$OUTDIR/top5_motif.txt

# Also extract top 10 for comparison
python3 src/extract_top_conserved.py \
  $SCRATCH/$OUTDIR/${PDB_ID_LOWER}_conservation.json \
  --top-n 10 \
  --exclude-gaps \
  --min-identity 0.2 \
  --mapping $SCRATCH/$OUTDIR/alignment_mapping.json \
  --output $SCRATCH/$OUTDIR/top10_motif.txt

#####################
# STEP 7: Benchmark against M-CSA (BASELINE)
#####################
echo ""
echo "[7] Benchmarking against M-CSA ground truth..."

if [ -f "$MCSA_FILE" ]; then
  # Build P2Rank argument if available
  P2RANK_ARG=""
  if [ -n "$P2RANK_JSON" ] && [ -f "$P2RANK_JSON" ]; then
    P2RANK_ARG="--p2rank-json $P2RANK_JSON"
  fi

  python3 src/benchmark_mcsa.py \
    $SCRATCH/$OUTDIR/${PDB_ID_LOWER}_conservation.json \
    $MCSA_FILE \
    $SCRATCH/$OUTDIR/alignment_mapping.json \
    --pdb-id $PDB_ID_LOWER \
    --pdb-file $PDB_CACHE/${PDB_ID}.pdb \
    --conservation-threshold $SCORE_T \
    --exclude-gaps \
    $P2RANK_ARG \
    --min-identity 0.2 \
    --output $SCRATCH/$OUTDIR/baseline_performance.json
  
  echo "  ✓ Baseline performance saved to: baseline_performance.json"
else
  echo "  WARNING: M-CSA file not found at $MCSA_FILE"
  echo "  Skipping benchmarking. Download from M-CSA database to enable."
fi

#####################
# STEP 8: Summary (updated)
#####################
echo ""
echo "[8] Done!"
echo "============================================"
echo "Output files in $SCRATCH/$OUTDIR/"
echo "  foldseek_hits.tsv         - raw foldseek easy-search output"
echo "  homologs.txt              - filtered UniProt IDs (foldseek hits, top-N)"
echo "  structures/               - downloaded AlphaFold PDBs"
echo "  foldmason_result_aa.fa    - structural MSA (amino acid)"
echo "  foldmason_result_3di.fa   - structural MSA (3Di)"
echo "  foldmason_result.nw       - guide tree"
echo "  ${PDB_ID_LOWER}_conservation.json  - conservation scores"
echo "  alignment_mapping.json    - alignment column → residue ID mapping"
echo "  top5_motif.txt            - top 5 conserved positions"
echo "  top10_motif.txt           - top 10 conserved positions"
echo "  baseline_performance.json - M-CSA benchmark results (if available)"
echo "============================================"


PIPELINE_END=$(date +%s)
PIPELINE_ELAPSED=$((PIPELINE_END - PIPELINE_START))
echo "Total time: ${PIPELINE_ELAPSED}s"
echo "ELAPSED_SECONDS=$PIPELINE_ELAPSED"