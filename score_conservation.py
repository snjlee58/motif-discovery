# python3 score_conservation.py $SCRATCH/foldmason_result_aa.fa "$SCRATCH/conservation_scores"

import sys
from pathlib import Path
from Bio import AlignIO
from evomotif_conservation import ConservationScorer

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

def read_fasta_alignment(path: str):
    headers, seqs = [], [] # list of names and sequences
    h = None # current sequence header '>'
    buf = [] # temporary storage for sequence lines (spanning multiple lines)
    with open(path) as f:
        for line in f:
            line = line.strip() # removes leading/trailing whitespace
            if not line:
                continue
            if line.startswith(">"): # join and store the previous sequence
                if h is not None:
                    seq = "".join(buf).upper()
                    headers.append(h)
                    seqs.append(seq)
                h = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if h is not None: # store the last sequence
            seq = "".join(buf).upper()
            headers.append(h)
            seqs.append(seq)

    if not seqs:
        raise ValueError(f"No sequences found in {path}.")

    L = len(seqs[0])
    for i, s in enumerate(seqs):
        if len(s) != L:
            raise ValueError(f"Not all sequences same length. Seq {i} has {len(s)} vs {L}")
    return headers, seqs

def main():
    outdir = Path(sys.argv[1])
    msa_path = Path(sys.argv[2]) if len(sys.argv) > 2 else outdir / "foldmason_result_aa.fa"
    protein_name = sys.argv[3] if len(sys.argv) > 3 else "1btl"
    cons_file = outdir / f"{protein_name}_conservation.json"

    # Step 1: Read the AA MSA
    print(f"\n[1] Reading alignment from {msa_path}...")
    alignment = AlignIO.read(msa_path, "fasta")
    print(f"       {len(alignment)} sequences, {alignment.get_alignment_length()} columns")

    scorer = ConservationScorer()

    # Step 2: Calculate AA conservation
    print(f"\n[2] Calculating conservation scores...")
    conservation = scorer.calculate_combined_conservation(alignment)
    print(f"       Mean: {conservation.mean():.3f}, Max: {conservation.max():.3f}")

    # Step 3: Calculate 3Di structural conservation (if available)
    di3_scores = None
    di3_path = outdir / "foldmason_result_3di.fa"
    if di3_path.exists():
        print(f"\n[3] Calculating 3Di structural conservation from {di3_path.name}...")
        di3_alignment = AlignIO.read(str(di3_path), "fasta")
        di3_scores = scorer.calculate_3di_conservation(di3_alignment)
        print(f"       Mean: {di3_scores.mean():.3f}, Max: {di3_scores.max():.3f}")
    else:
        print(f"\n[3] 3Di MSA not found at {di3_path}, skipping...")

    scorer.save_conservation_scores(alignment, cons_file, di3_conservation=di3_scores)
    print(f"\n       Saved to {cons_file}")


if __name__ == "__main__":
    main()