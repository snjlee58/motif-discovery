#!/bin/bash
set -euo pipefail

# Pipeline: query → find AFDB cluster → get all members → foldmason

# Map query protein to afdb cluster representative

# Get pdb of members of that cluster

# Foldmason MSA of the members

# Conservation scoring

# for each protein:
#     foldseek easy-search -> foldmason

# for each cluster:
#     foldmason easy-msa cluster_members/
    # foldmason easy-msa beta_lactamase_cluster/ msaDB msa_result tmp

# Prerequisites:
#   1. Download AFDB cluster membership file:
#  wget https://afdb-cluster.steineggerlab.workers.dev/1-AFDBClusters-entryId_repId_taxId.tsv.gz \
#    -O $SCRATCH/afdb_clusters/1-AFDBClusters-entryId_repId_taxId.tsv.gz
#  gunzip $SCRATCH/afdb_clusters/1-AFDBClusters-entryId_repId_taxId.tsv.gz
#   2. (Optional) If you want ALL members including non-AFDB50:
# wget https://afdb-cluster.steineggerlab.workers.dev/5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz \
#   -O $SCRATCH/afdb_clusters/5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz
# gunzip $SCRATCH/afdb_clusters/5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz

# Family View Pipeline
# Instead of: query → foldseek search → hits → foldmason
# We do:      query → find AFDB cluster → get all members → foldmason
#
# Prerequisites:
#   1. Download AFDB cluster membership file:
#      wget https://afdb-cluster.steineggerlab.workers.dev/5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz \
#        -O $FAST/afdb_clusters/5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz
#      gunzip $FAST/afdb_clusters/5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz
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

PDB_ID_LOWER=$(echo "$PDB_ID" | tr '[:upper:]' '[:lower:]')

# $FAST points to shared NVMe storage (e.g. /fast/sunny) holding persistent reference
# databases (AFDB clusters, M-CSA, PDB cache) — visible from every compute node.
# Falls back to $SCRATCH for single-node setups where nothing has been moved yet.
FAST="${FAST:-$SCRATCH}"

CLUSTER_FILE="$FAST/afdb_clusters/5-allmembers-repId-entryId-cluFlag-taxId.tsv"
MCSA_FILE="$FAST/m-csa/catalytic_residues_homologues_parsed.tsv"
PDB_CACHE="$FAST/pdb_files"
mkdir -p "$PDB_CACHE"

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
echo "Family View Pipeline"
echo "  UniProt:    $UNIPROT_ID"
echo "  PDB:        $PDB_ID"
echo "  Output:     $SCRATCH/$OUTDIR"
echo "============================================"

#####################
# STEP 1: Find AFDB cluster for this protein
#####################
echo ""
echo "[1] Finding AFDB cluster for $UNIPROT_ID..."

if [ ! -f "$CLUSTER_FILE" ]; then
  echo "ERROR: Cluster file not found at $CLUSTER_FILE"
  echo "Download it first (FAST=$FAST):"
  echo "  mkdir -p $FAST/afdb_clusters"
  echo "  wget https://afdb-cluster.steineggerlab.workers.dev/5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz -O $FAST/afdb_clusters/5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz"
  echo "  gunzip $FAST/afdb_clusters/5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz"
  exit 1
fi

# # File 1 format: entryId \t repId \t taxId
# # entryId is a UniProt accession
# # Find the cluster representative for our protein
# REP_ID=$(grep -m1 "^${UNIPROT_ID}\b" "$CLUSTER_FILE" | cut -f2)

# if [ -z "$REP_ID" ]; then
#   echo "ERROR: $UNIPROT_ID not found in AFDB clusters"
#   echo "  Try searching with a different UniProt ID"
#   exit 1
# fi

# echo "  Cluster representative: $REP_ID"

# # Get all members in this cluster (all entries with the same repId)
# grep "\b${REP_ID}\b" "$CLUSTER_FILE" | cut -f1 > $SCRATCH/$OUTDIR/cluster_members_all.txt
# N_TOTAL=$(wc -l < $SCRATCH/$OUTDIR/cluster_members_all.txt)
# echo "  Total cluster members: $N_TOTAL"



# File 5 format: repId \t memId \t cluFlag \t taxId
# Find the cluster representative for our protein (our protein is in col 2)
REP_ID=$(grep "${UNIPROT_ID}" "$CLUSTER_FILE" | head -1 | cut -f1) || true

if [ -z "$REP_ID" ]; then
  echo "ERROR: $UNIPROT_ID not found in AFDB clusters"
  echo "  Try searching with a different UniProt ID"
  exit 1
fi

echo "  Cluster representative: $REP_ID"

# Get all members of that cluster (representative is in col 1)
grep "^${REP_ID}" "$CLUSTER_FILE" | cut -f2 > $SCRATCH/$OUTDIR/cluster_members_all.txt
N_TOTAL=$(wc -l < $SCRATCH/$OUTDIR/cluster_members_all.txt)
echo "  Total cluster members: $N_TOTAL"

cp $SCRATCH/$OUTDIR/cluster_members_all.txt $SCRATCH/$OUTDIR/cluster_members.txt
echo "  Using $N_TOTAL members"

#####################
# STEP 2: Download AlphaFold structures for cluster members
#####################
echo ""
echo "[2] Downloading AlphaFold structures..."
mkdir -p $SCRATCH/$OUTDIR/structures

TOTAL=$(wc -l < $SCRATCH/$OUTDIR/cluster_members.txt)

# Build the to-do list (skip already-cached structures from prior runs).
TODO_FILE="$SCRATCH/$OUTDIR/to_download.txt"
> "$TODO_FILE"
while read -r MEMBER_ID; do
  OUTFILE="$SCRATCH/$OUTDIR/structures/AF-${MEMBER_ID}-F1-model_v6.pdb"
  [ -f "$OUTFILE" ] || echo "$MEMBER_ID" >> "$TODO_FILE"
done < $SCRATCH/$OUTDIR/cluster_members.txt
NEED=$(wc -l < "$TODO_FILE")
CACHED=$((TOTAL - NEED))
echo "  $TOTAL total, $CACHED already cached, $NEED to download"

# Parallel downloads. AFDB EBI handles tens of concurrent connections fine,
# but if the surrounding batch (run_full_pipeline.sh) is also parallelizing
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

# Add the experimental PDB structure for the query
PDB_FILE="$PDB_CACHE/${PDB_ID}.pdb"
if [ ! -f "$PDB_FILE" ]; then
  echo "  Downloading $PDB_ID from RCSB..."
  wget -q --timeout=30 \
    "https://files.rcsb.org/download/${PDB_ID}.pdb" \
    -O "$PDB_FILE" 2>/dev/null || {
    echo "  WARNING: Could not download $PDB_ID.pdb from RCSB"
  }
fi

if [ -f "$PDB_FILE" ]; then
  echo "  Copying $PDB_ID.pdb into structures/"
  cp "$PDB_FILE" $SCRATCH/$OUTDIR/structures/
else
  echo "  WARNING: No experimental PDB for $PDB_ID. Will use AlphaFold structure in MSA."
fi

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
python3 score_conservation.py $SCRATCH/$OUTDIR $MSA_FILE $PDB_ID_LOWER

#####################
# STEP 5: Create alignment mapping
#####################
echo ""
echo "[5] Mapping alignment columns to PDB residue IDs..."
python3 map_alignment_to_pdb.py \
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
  if ! python3 parse_p2rank.py "$csv" -o "$P2RANK_JSON"; then
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
python3 extract_top_conserved.py \
  $SCRATCH/$OUTDIR/${PDB_ID_LOWER}_conservation.json \
  --top-n 5 \
  --exclude-gaps \
  --min-identity 0.2 \
  --mapping $SCRATCH/$OUTDIR/alignment_mapping.json \
  --output $SCRATCH/$OUTDIR/top5_motif.txt

# Also extract top 10 for comparison
python3 extract_top_conserved.py \
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

  python3 benchmark_mcsa.py \
    $SCRATCH/$OUTDIR/${PDB_ID_LOWER}_conservation.json \
    $MCSA_FILE \
    $SCRATCH/$OUTDIR/alignment_mapping.json \
    --pdb-id $PDB_ID_LOWER \
    --pdb-file $PDB_CACHE/${PDB_ID}.pdb \
    --top-n auto \
    --exclude-gaps \
    --catalytic-propensity \
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
echo "  cluster_members.txt       - UniProt IDs in this cluster"
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