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


    # Step 1: Read the MSA
    print(f"\n[1] 📊 Reading alignment from {msa_path}...")
    alignment = AlignIO.read(msa_path, "fasta")  # Converts Foldmason MSA fasta into MultipleSeqAlignment object
    # print(alignment)
    print(f"       {len(alignment)} sequences, {alignment.get_alignment_length()} columns")


    # Step 2: Calculate conservation
    print(f"\n[2] 📊 Calculating conservation scores...")
    scorer = ConservationScorer()
    
    conservation = scorer.calculate_combined_conservation(alignment)
    gap_freq = scorer.calculate_gap_frequency(alignment)
    consensus = scorer.get_consensus_sequence(alignment)
    
    scorer.save_conservation_scores(alignment, cons_file)
    print(f"       Mean conservation: {conservation.mean():.3f}, "
          f"Max: {conservation.max():.3f}")
    print(f"       Saved to {cons_file}")
    
    # results_data['conservation'] = conservation.tolist()
    # results_data['mean_conservation'] = float(conservation.mean())
    # results_data['max_conservation'] = float(conservation.max())
    # results_data['files']['conservation'] = cons_file


if __name__ == "__main__":
    main()