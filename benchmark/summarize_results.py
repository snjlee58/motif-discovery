#!/usr/bin/env python3
"""
Summarize benchmark results across all proteins.

Usage (from anywhere):
    python3 ~/motif/benchmark/summarize_results.py \
        ~/motif/benchmark/mcsa_representatives_parsed_monomers.tsv \
        --scratch /fast/sunny/motif/batch_runs/<TS>_job<ID>
"""
import json
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('tsv', help='Protein list TSV (mcsa_id, pdb_id, n_catalytic_residues)')
    parser.add_argument('--scratch', default='/mnt/scratch/sunny',
                        help='Directory containing per-PDB result subdirs')
    parser.add_argument('--result-file', default='baseline_performance.json',
                    help='Which result file to read (e.g. performance_top10.json)')
    args = parser.parse_args()

    scratch = Path(args.scratch)

    proteins = []
    with open(args.tsv) as f:
        next(f)  # skip header
        for line in f:
            mcsa_id, pdb, n_res = line.strip().split('\t')
            proteins.append({
                'mcsa_id': mcsa_id,
                'pdb': pdb.upper(),
                'n_res': int(n_res)
            })

    print(f"\n{'PDB':<8} {'MCSA':<8} {'Cat.Res':<9} "
          f"{'Precision':<11} {'Recall':<9} {'F1':<8} "
          f"{'TP/True':<10} {'Top-N':<7} {'Status'}")
    print("=" * 86)

    results = []
    not_run = []

    for p in proteins:
        pdb = p['pdb']
        matches = sorted(scratch.glob(f"{pdb}/{args.result_file}"))

        if not matches:
            print(f"{pdb:<8} {p['mcsa_id']:<8} {p['n_res']:<9} {'NOT RUN'}")
            not_run.append(pdb)
            continue

        with open(matches[-1]) as f:
            r = json.load(f)

        m = r['metrics']
        results.append(m)

        print(f"{pdb:<8} {p['mcsa_id']:<8} {p['n_res']:<9} "
              f"{m['precision']:<11.3f} {m['recall']:<9.3f} {m['f1']:<8.3f} "
              f"{m['tp']}/{m['n_true']:<8} {r['top_n']:<7} ✓")

    if results:
        print("=" * 86)
        n = len(results)

        def avg(key):
            return sum(r[key] for r in results) / n

        print(f"\nSummary ({n}/{len(proteins)} proteins completed):")
        print(f"  Mean Precision: {avg('precision'):.3f}")
        print(f"  Mean Recall:    {avg('recall'):.3f}")
        print(f"  Mean F1:        {avg('f1'):.3f}")
        print(f"  Min F1:         {min(r['f1'] for r in results):.3f}")
        print(f"  Max F1:         {max(r['f1'] for r in results):.3f}")

    if not_run:
        print(f"\nNot yet run: {not_run}")

if __name__ == '__main__':
    main()