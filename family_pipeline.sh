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

# Usage:
#   bash family_pipeline.sh <pdb_id> [uniprot_id] [output_dir]
#   If uniprot_id is omitted, it will be resolved automatically via UniProt ID Mapping API.
#
# Examples:
#   bash family_pipeline.sh 1BTL
#   bash family_pipeline.sh 1BTL P62593
#   bash family_pipeline.sh 1BTL P62593 250212_family_1btl



# Family View Pipeline
# Instead of: query → foldseek search → hits → foldmason
# We do:      query → find AFDB cluster → get all members → foldmason
#
# Prerequisites:
#   1. Download AFDB cluster membership file:
#      wget https://afdb-cluster.steineggerlab.workers.dev/1-AFDBClusters-entryId_repId_taxId.tsv.gz \
#        -O $SCRATCH/afdb_clusters/1-AFDBClusters-entryId_repId_taxId.tsv.gz
#      gunzip $SCRATCH/afdb_clusters/1-AFDBClusters-entryId_repId_taxId.tsv.gz
#
#   2. (Optional) If you want ALL members including non-AFDB50:
#      wget https://afdb-cluster.steineggerlab.workers.dev/5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz
#
# Usage:
#   bash family_pipeline.sh <pdb_id> [uniprot_id] [output_dir]
#
# Examples:
#   bash family_pipeline.sh 1BTL
#   bash family_pipeline.sh 1BTL P62593
#   bash family_pipeline.sh 1BTL P62593 250212_family_1btl

#####################
# ARGUMENTS
#####################
PDB_ID=${1:?"Usage: bash family_pipeline.sh <pdb_id> [uniprot_id] [output_dir]"}
UNIPROT_ID=${2:-""}
OUTDIR=${3:-$(date +%y%m%d_%H%M%S)_family_${PDB_ID}}

PDB_ID_LOWER=$(echo "$PDB_ID" | tr '[:upper:]' '[:lower:]')
# CLUSTER_FILE="$SCRATCH/afdb_clusters/1-AFDBClusters-entryId_repId_taxId.tsv"
CLUSTER_FILE="$SCRATCH/afdb_clusters/5-allmembers-repId-entryId-cluFlag-taxId.tsv"
MAX_MEMBERS=100  # cap to keep FoldMason tractable
P2RANK_DIR="$SCRATCH/p2rank_2.5"  # Path to P2Rank installation

#####################
# STEP 0: Resolve UniProt ID if not provided
#####################
if [ -z "$UNIPROT_ID" ]; then
  echo "Resolving UniProt ID for PDB: $PDB_ID via UniProt ID Mapping API..."

  # Step 1: Submit mapping job (PDB -> UniProtKB)
  SUBMIT_RESPONSE=$(curl -s --form 'from="PDB"' --form 'to="UniProtKB"' \
    --form "ids=\"${PDB_ID_LOWER}\"" \
    https://rest.uniprot.org/idmapping/run)
  echo "  Submit response: $SUBMIT_RESPONSE"

  JOB_ID=$(echo "$SUBMIT_RESPONSE" | python3 -c "import sys,json; print(json.load(sys.stdin)['jobId'])" 2>/dev/null)

  if [ -z "$JOB_ID" ]; then
    echo "ERROR: Failed to submit UniProt ID mapping job"
    exit 1
  fi
  echo "  Submitted job: $JOB_ID"

  # Step 2: Poll until complete, following redirects (-L)
  for i in $(seq 1 15); do
    sleep 1
    RESULT_RAW=$(curl -s -L "https://rest.uniprot.org/idmapping/status/${JOB_ID}")

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
    UNIPROT_ID=$(curl -s -L "https://rest.uniprot.org/idmapping/results/stream/${JOB_ID}" \
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
    echo "  Try providing it manually: bash family_pipeline.sh $PDB_ID <uniprot_id>"
    exit 1
  fi
  echo "  Resolved: $PDB_ID -> $UNIPROT_ID"
fi

mkdir -p $SCRATCH/$OUTDIR
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
  echo "Download it first:"
  echo "  mkdir -p $SCRATCH/afdb_clusters"
  echo "  wget https://afdb-cluster.steineggerlab.workers.dev/1-AFDBClusters-entryId_repId_taxId.tsv.gz -O $SCRATCH/afdb_clusters/1-AFDBClusters-entryId_repId_taxId.tsv.gz"
  echo "  gunzip $SCRATCH/afdb_clusters/1-AFDBClusters-entryId_repId_taxId.tsv.gz"
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




# Cap the number of members
if [ "$N_TOTAL" -gt "$MAX_MEMBERS" ]; then
  echo "  Capping to $MAX_MEMBERS members (random sample + ensuring query is included)"
  # Keep query, randomly sample the rest
  grep -v "^${UNIPROT_ID}$" $SCRATCH/$OUTDIR/cluster_members_all.txt | shuf | head -n $((MAX_MEMBERS - 1)) > $SCRATCH/$OUTDIR/cluster_members.txt || true
  echo "$UNIPROT_ID" >> $SCRATCH/$OUTDIR/cluster_members.txt
else
  cp $SCRATCH/$OUTDIR/cluster_members_all.txt $SCRATCH/$OUTDIR/cluster_members.txt
fi

N_MEMBERS=$(wc -l < $SCRATCH/$OUTDIR/cluster_members.txt)
echo "  Using $N_MEMBERS members"

#####################
# STEP 2: Download AlphaFold structures for cluster members
#####################
echo ""
echo "[2] Downloading AlphaFold structures..."
mkdir -p $SCRATCH/$OUTDIR/structures

DOWNLOADED=0
FAILED=0
FAILED_IDS=""

while read -r MEMBER_ID; do
  OUTFILE="$SCRATCH/$OUTDIR/structures/AF-${MEMBER_ID}-F1-model_v6.pdb"
  
  if [ -f "$OUTFILE" ]; then
    DOWNLOADED=$((DOWNLOADED + 1))
    continue
  fi

  URL="https://alphafold.ebi.ac.uk/files/AF-${MEMBER_ID}-F1-model_v6.pdb"
  
  if wget -q --timeout=10 -O "$OUTFILE" "$URL" 2>/dev/null; then
    DOWNLOADED=$((DOWNLOADED + 1))
  else
    rm -f "$OUTFILE"
    FAILED=$((FAILED + 1))
    FAILED_IDS="${FAILED_IDS} ${MEMBER_ID}"
  fi
done < $SCRATCH/$OUTDIR/cluster_members.txt

echo "  Downloaded: $DOWNLOADED"
echo "  Failed: $FAILED"
if [ -n "$FAILED_IDS" ]; then
  echo "  Failed IDs:$FAILED_IDS"
fi

echo "  DEBUG: finished download loop"

# Add the experimental PDB structure for the query
PDB_FILE="$SCRATCH/${PDB_ID}.pdb"
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

# # Add the experimental PDB structure for the query
# cp $SCRATCH/${PDB_ID}.pdb $SCRATCH/$OUTDIR/structures/ 2>/dev/null || true

# N_STRUCTURES=$(ls $SCRATCH/$OUTDIR/structures/*.pdb $SCRATCH/$OUTDIR/structures/*.cif 2>/dev/null | wc -l)
# echo "  Total structures for alignment: $N_STRUCTURES"

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
  $SCRATCH/tmp \
  --report-mode 1

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
python3 pipeline.py $SCRATCH/$OUTDIR $MSA_FILE $PDB_ID_LOWER

#####################
# STEP 5: Create alignment mapping
#####################
echo ""
echo "[5] Mapping alignment columns to PDB residue IDs..."
python3 map_alignment_to_pdb.py \
  $MSA_FILE \
  $PDB_ID \
  --pdb-file $SCRATCH/${PDB_ID}.pdb \
  --uniprot $UNIPROT_ID \
  -o $SCRATCH/$OUTDIR/alignment_mapping.json

#####################
# STEP 5b: P2Rank binding site prediction
#####################
echo ""
echo "[5b] Running P2Rank binding site prediction..."

P2RANK_JSON="$SCRATCH/$OUTDIR/p2rank_scores.json"
PDB_FILE="$SCRATCH/${PDB_ID}.pdb"

if [ -d "$P2RANK_DIR" ] && [ -f "$PDB_FILE" ]; then
  $P2RANK_DIR/prank predict -f "$PDB_FILE" -o $SCRATCH/$OUTDIR/p2rank_output 2>/dev/null

  # Find the residues CSV
  P2RANK_CSV=$(find $SCRATCH/$OUTDIR/p2rank_output -name "*_residues.csv" | head -1)

  if [ -n "$P2RANK_CSV" ]; then
    python3 parse_p2rank.py "$P2RANK_CSV" -o "$P2RANK_JSON"
    echo "  ✓ P2Rank scores saved to: p2rank_scores.json"
  else
    echo "  WARNING: P2Rank ran but no residues CSV found"
    P2RANK_JSON=""
  fi
else
  echo "  WARNING: P2Rank not found at $P2RANK_DIR or PDB file missing"
  echo "  Skipping P2Rank. Install: cd \$SCRATCH && wget https://github.com/rdk/p2rank/releases/download/2.5/p2rank_2.5.tar.gz && tar xzf p2rank_2.5.tar.gz"
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

MCSA_FILE="$SCRATCH/m-csa/catalytic_residues_separated_aligned.tsv"

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