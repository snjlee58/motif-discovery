#!/usr/bin/env python3
# diagnose_benchmark.py
# Usage: python3 diagnose_benchmark.py selected_proteins.tsv
import json, argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('tsv')
    parser.add_argument('--scratch', default='/mnt/scratch/sunny')
    args = parser.parse_args()
    scratch = Path(args.scratch)

    proteins = []
    with open(args.tsv) as f:
        next(f)
        for line in f:
            mcsa_id, uniprot, pdb, n_res = line.strip().split('\t')
            proteins.append({'pdb': pdb.upper(), 'uniprot': uniprot})

    for p in proteins:
        pdb = p['pdb']
        pdb_lower = pdb.lower()

        # Find output dir
        perf_matches = sorted(scratch.glob(f"*family_{pdb}*/baseline_performance.json"))
        cons_matches = sorted(scratch.glob(f"*family_{pdb}*/{pdb_lower}_conservation.json"))
        map_matches  = sorted(scratch.glob(f"*family_{pdb}*/alignment_mapping.json"))

        if not perf_matches:
            print(f"\n{pdb}: no baseline_performance.json found")
            continue

        with open(perf_matches[-1]) as f:
            perf = json.load(f)
        with open(cons_matches[-1]) as f:
            cons = json.load(f)
        with open(map_matches[-1]) as f:
            mapping = json.load(f)

        m = perf['metrics']
        true_set = set(perf['mcsa_ground_truth'])
        predicted_set = set(perf['predicted'])
        n_seqs = cons['n_sequences']
        aln_len = cons['alignment_length']

        # Build resid -> conservation rank lookup
        # mapping: alignment_col(str) -> resid
        resid_to_col = {v: int(k) for k, v in mapping['mapping'].items()}
        col_to_pos = {p['position']: p for p in cons['positions']}

        print(f"\n{'='*60}")
        print(f"{pdb} | F1={m['f1']:.3f} | {n_seqs} seqs | aln_len={aln_len}")
        print(f"  Ground truth ({len(true_set)}): {sorted(true_set)}")
        print(f"  Predicted:    {sorted(predicted_set)}")
        print(f"\n  Catalytic residue ranks (out of {aln_len} columns):")

        for resid in sorted(true_set):
            col = resid_to_col.get(resid)
            if col is None:
                print(f"    ResID {resid:>4}: NOT IN MAPPING (gap in query?)")
                continue
            pos = col_to_pos.get(col)
            if pos is None:
                print(f"    ResID {resid:>4}: col {col} NOT IN CONSERVATION DATA")
                continue

            # Get rank
            all_cons = sorted(cons['positions'],
                              key=lambda x: x['conservation'], reverse=True)
            rank = next((i+1 for i, p in enumerate(all_cons)
                        if p['position'] == col), -1)
            hit = "✓" if resid in predicted_set else "✗"
            print(f"    {hit} ResID {resid:>4} ({pos['consensus']}) | "
                  f"col={col:>4} | cons={pos['conservation']:.3f} | "
                  f"identity={pos['identity']:.3f} | rank={rank}/{len(all_cons)}")

if __name__ == '__main__':
    main()