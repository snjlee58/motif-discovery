#!/usr/bin/env python3
"""
Diagnose a benchmark TSV: separate ranking failures from pool-dropout failures,
print F1 distribution, and find ranking-failure cases worth investigating.

Usage:
    python3 analysis/diagnose_benchmark.py <benchmark_summary.tsv> [batch_dir]

If batch_dir is given, also reads per-PDB JSONs to find ranking failures
(catalytic residues that ARE in the conservation pool but ranked below top-K).
"""
import csv, json, statistics, sys
from pathlib import Path


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    tsv = Path(sys.argv[1])
    batch_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else None

    rows = []
    with open(tsv) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            rows.append({
                'pdb': r['pdb_id'],
                'f1': float(r['f1']),
                'n_pred': int(r['n_predicted']),
                'n_true': int(r['n_true']),
            })

    n = len(rows)
    f1s = [r['f1'] for r in rows]

    print(f"\n=== F1 distribution (n={n}) ===")
    print(f"  mean   {statistics.mean(f1s):.3f}")
    print(f"  median {statistics.median(f1s):.3f}")
    print(f"  std    {statistics.stdev(f1s):.3f}")
    buckets = [(0.0, 0.0001), (0.0001, 0.2), (0.2, 0.4), (0.4, 0.6),
               (0.6, 0.8), (0.8, 1.0001)]
    for lo, hi in buckets:
        cnt = sum(1 for x in f1s if lo <= x < hi)
        label = f"{lo:.2f}<=F1<{hi:.2f}" if lo > 0 else "F1==0"
        print(f"  {label:<16} {cnt:4d}  ({100*cnt/n:.1f}%)")
    print(f"  F1==1.0          {sum(1 for x in f1s if x == 1.0):4d}")

    # === Pool dropout vs ranking failure ===
    dropout_total = [r for r in rows if r['n_pred'] == 0]
    dropout_partial = [r for r in rows if 0 < r['n_pred'] < r['n_true']]
    full_pred = [r for r in rows if r['n_pred'] == r['n_true']]

    print(f"\n=== Failure mode breakdown ===")
    print(f"  Total dropout   (n_pred=0):              {len(dropout_total):4d}  "
          f"-> filter killed everything; F1 mean = "
          f"{statistics.mean(r['f1'] for r in dropout_total) if dropout_total else 0:.3f}")
    print(f"  Partial dropout (0 < n_pred < n_true):    {len(dropout_partial):4d}  "
          f"-> filter killed some; F1 mean = "
          f"{statistics.mean(r['f1'] for r in dropout_partial) if dropout_partial else 0:.3f}")
    print(f"  Full pool       (n_pred = n_true):        {len(full_pred):4d}  "
          f"-> filter kept candidates; F1 mean = "
          f"{statistics.mean(r['f1'] for r in full_pred) if full_pred else 0:.3f}")

    if dropout_total:
        print(f"\n  Total-dropout PDBs: {[r['pdb'] for r in dropout_total]}")
    if dropout_partial:
        print(f"  Partial-dropout PDBs: {[(r['pdb'], r['n_pred'], r['n_true']) for r in dropout_partial]}")

    print(f"\n=== Performance among proteins where filter worked ===")
    if full_pred:
        ff = [r['f1'] for r in full_pred]
        print(f"  Mean F1 (filter-OK only):   {statistics.mean(ff):.3f}")
        print(f"  Median F1 (filter-OK only): {statistics.median(ff):.3f}")
        print(f"  F1==0 (filter-OK only):     {sum(1 for x in ff if x == 0):d} "
              f"({100*sum(1 for x in ff if x == 0)/len(ff):.1f}%)  <- these are PURE ranking failures")

    # === If batch_dir given, classify ranking failures: in-pool-but-low-rank vs not-in-conservation-data ===
    if not batch_dir or not batch_dir.exists():
        print("\n[skip] No batch_dir given; skipping rank-of-true-residue analysis.")
        return

    print(f"\n=== Rank-of-true-residue analysis (where F1==0 and n_pred==n_true) ===")
    pure_rank_fail = [r for r in full_pred if r['f1'] == 0]
    print(f"Examining {len(pure_rank_fail)} pure ranking failures...\n")

    in_pool_summary = []
    for r in pure_rank_fail[:30]:  # cap output
        pdb = r['pdb'].lower()
        cons_path = batch_dir / r['pdb'] / f"{pdb}_conservation.json"
        map_path  = batch_dir / r['pdb'] / "alignment_mapping.json"
        result_path = batch_dir / r['pdb'] / "baseline_performance.json"
        if not (cons_path.exists() and map_path.exists() and result_path.exists()):
            continue
        with open(cons_path) as f: cons = json.load(f)
        with open(map_path)  as f: mp   = json.load(f)['mapping']
        with open(result_path) as f: res = json.load(f)
        true_set = set(res['mcsa_ground_truth'])

        # Build a sorted-by-conservation list of (resid, conservation)
        col_to_resid = {int(k): v for k, v in mp.items()}
        scored = []
        for p in cons['positions']:
            resid = col_to_resid.get(p['position'])
            if resid is None: continue
            scored.append((resid, p['conservation']))
        scored.sort(key=lambda x: x[1], reverse=True)
        ranks = {resid: i+1 for i, (resid, _) in enumerate(scored)}

        true_ranks = [(t, ranks.get(t)) for t in true_set]
        in_pool = [(t, rk) for t, rk in true_ranks if rk is not None]
        not_in  = [t for t, rk in true_ranks if rk is None]
        n_total = len(scored)

        print(f"  {r['pdb']}: K={r['n_true']}, pool={n_total}")
        for t, rk in sorted(true_ranks, key=lambda x: (x[1] is None, x[1])):
            if rk is None:
                print(f"     - res {t}: NOT in conservation pool (filter dropped, or chain mismatch)")
            else:
                pct = 100 * rk / n_total
                print(f"     - res {t}: rank {rk}/{n_total} ({pct:.1f}%)")
        in_pool_summary.append({'pdb': r['pdb'], 'in_pool': len(in_pool), 'not_in': len(not_in),
                                'n_true': r['n_true']})

    if in_pool_summary:
        print(f"\n--- Summary for pure-ranking-failure cases examined ---")
        total_true = sum(s['n_true'] for s in in_pool_summary)
        total_in_pool = sum(s['in_pool'] for s in in_pool_summary)
        total_not_in = sum(s['not_in'] for s in in_pool_summary)
        print(f"  Across {len(in_pool_summary)} PDBs, {total_true} true catalytic residues:")
        print(f"     {total_in_pool} ({100*total_in_pool/total_true:.0f}%) in conservation pool but ranked too low")
        print(f"     {total_not_in} ({100*total_not_in/total_true:.0f}%) not in conservation pool at all")


if __name__ == '__main__':
    main()
