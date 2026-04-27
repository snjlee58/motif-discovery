#!/usr/bin/env python3
"""
Compare benchmark performance across multiple top-N values.

Reads the per-PDB JSONs produced by `rescore_batch.sh` with TOP_N_MULTIPLIERS
set, and prints a side-by-side comparison: precision / recall / F1 at each
multiplier, plus aggregate stats.

Lets you separate two failure modes:
  - "missing the residues entirely" → recall stays flat as N grows
  - "ranking them too low"          → recall improves with N, precision drops

Usage:
    python3 analysis/analyze_top_n_comparison.py <rescore_dir>

Example:
    python3 analysis/analyze_top_n_comparison.py \
        /fast/sunny/motif/batch_runs/260426_025725_job527128/rescore_260427_113000_topn1-2-3
"""

import argparse
import json
import sys
from pathlib import Path
from statistics import mean, median


def load_results(rescore_dir):
    """Yield (pdb_id, runs_dict, n_true) for each PDB JSON in the dir.

    runs_dict: {label: {top_n, predicted, metrics}}.
    Falls back to single-value 'primary' when the JSON has no metrics_by_top_n.
    """
    for jf in sorted(Path(rescore_dir).glob('*.json')):
        if jf.name in ('benchmark_summary.json',):
            continue
        try:
            with open(jf) as f:
                d = json.load(f)
        except json.JSONDecodeError:
            continue
        if 'metrics_by_top_n' in d:
            yield d['pdb_id'], d['metrics_by_top_n'], d.get('n_true')
        elif 'metrics' in d:
            yield d['pdb_id'], {'primary': {'top_n': d.get('top_n'),
                                            'predicted': d.get('predicted', []),
                                            'metrics': d['metrics']}}, None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('rescore_dir', help='Rescore output directory '
                        '(<batch_dir>/rescore_<TS>_topn<...>/)')
    parser.add_argument('--per-pdb', action='store_true',
                        help='Print the per-PDB row table (default: just aggregate)')
    args = parser.parse_args()

    rescore_dir = Path(args.rescore_dir)
    if not rescore_dir.is_dir():
        sys.exit(f"ERROR: not a directory: {rescore_dir}")

    rows = list(load_results(rescore_dir))
    if not rows:
        sys.exit(f"ERROR: no per-PDB JSONs found in {rescore_dir}")

    # Discover labels (preserve insertion order from the first row)
    labels = list(rows[0][1].keys())
    n_pdbs = len(rows)

    print(f"Loaded {n_pdbs} PDB results from {rescore_dir}")
    print(f"Top-N labels compared: {labels}\n")

    # Per-PDB table (optional)
    if args.per_pdb:
        head = f"{'PDB':<8}{'n_true':>7}  "
        for lab in labels:
            head += f"{'P@'+lab:>7}{'R@'+lab:>7}{'F1@'+lab:>7}  "
        print(head)
        print("=" * len(head))
        for pdb, runs, n_true in rows:
            line = f"{pdb:<8}{(n_true or '-'):>7}  "
            for lab in labels:
                m = runs[lab]['metrics']
                line += f"{m['precision']:>7.3f}{m['recall']:>7.3f}{m['f1']:>7.3f}  "
            print(line)
        print()

    # Aggregate stats per label
    print(f"Aggregate stats across {n_pdbs} PDBs")
    print("=" * 70)
    print(f"  {'label':<8}{'mean P':>9}{'mean R':>9}{'mean F1':>10}"
          f"{'med P':>9}{'med R':>9}{'med F1':>9}{'F1=0':>7}{'F1≥.5':>7}")
    for lab in labels:
        ps, rs, f1s = [], [], []
        for _, runs, _ in rows:
            m = runs[lab]['metrics']
            ps.append(m['precision']); rs.append(m['recall']); f1s.append(m['f1'])
        zeros = sum(1 for v in f1s if v == 0)
        good = sum(1 for v in f1s if v >= 0.5)
        print(f"  {lab:<8}{mean(ps):>9.3f}{mean(rs):>9.3f}{mean(f1s):>10.3f}"
              f"{median(ps):>9.3f}{median(rs):>9.3f}{median(f1s):>9.3f}"
              f"{zeros:>7}{good:>7}")

    # Diagnostic: how many PDBs improve recall as N grows?
    if len(labels) >= 2:
        print()
        first, last = labels[0], labels[-1]
        improved, flat, worsened = 0, 0, 0
        for _, runs, _ in rows:
            r0 = runs[first]['metrics']['recall']
            r1 = runs[last]['metrics']['recall']
            if r1 > r0 + 0.001:
                improved += 1
            elif r1 < r0 - 0.001:
                worsened += 1
            else:
                flat += 1
        print(f"Recall change from {first} → {last}:")
        print(f"  improved : {improved:>4}  ({100*improved/n_pdbs:.1f}%)  "
              f"— catalytic residues are present, just ranked low")
        print(f"  flat     : {flat:>4}  ({100*flat/n_pdbs:.1f}%)  "
              f"— either already perfect at small N, or genuinely missing")
        print(f"  worsened : {worsened:>4}  ({100*worsened/n_pdbs:.1f}%)  "
              f"— shouldn't happen unless ties; investigate if non-zero")


if __name__ == '__main__':
    main()
