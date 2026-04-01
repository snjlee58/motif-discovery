#!/usr/bin/env python3
"""
Aggregate benchmark results from batch pipeline runs.

Reads batch_summary.tsv and produces overall statistics.

Usage:
    python3 aggregate_benchmark.py <batch_summary.tsv>
"""

import csv
import argparse
import json
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(description='Aggregate benchmark results')
    parser.add_argument('summary_tsv', help='batch_summary.tsv from batch_pipeline.sh')
    parser.add_argument('-o', '--output', default=None, help='Save detailed results as JSON')
    parser.add_argument('--top-n', type=int, nargs='+', default=[5, 10],
                        help='Top-N values to report (default: 5 10)')
    
    args = parser.parse_args()
    
    # Read results
    entries = []
    with open(args.summary_tsv, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            entries.append(row)
    
    total = len(entries)
    succeeded = [e for e in entries if e['status'] == 'SUCCESS']
    failed = [e for e in entries if e['status'] != 'SUCCESS']
    
    print("=" * 70)
    print("BENCHMARK SUMMARY")
    print("=" * 70)
    
    print(f"\nTotal entries:    {total}")
    print(f"Succeeded:        {len(succeeded)}")
    print(f"Failed:           {len(failed)}")
    if total > 0:
        print(f"Success rate:     {len(succeeded)/total*100:.1f}%")
    
    if not succeeded:
        print("\nNo successful runs to analyze.")
        return
    
    # Extract metrics
    precisions = []
    recalls = []
    f1s = []
    
    for e in succeeded:
        try:
            precisions.append(float(e['precision']))
            recalls.append(float(e['recall']))
            f1s.append(float(e['f1']))
        except (ValueError, KeyError):
            continue
    
    n = len(precisions)
    
    print(f"\n--- Performance Metrics (n={n}) ---")
    print(f"{'Metric':<15} {'Mean':>8} {'Median':>8} {'Std':>8} {'Min':>8} {'Max':>8}")
    print("-" * 60)
    
    for name, values in [('Precision', precisions), ('Recall', recalls), ('F1', f1s)]:
        if not values:
            continue
        mean = sum(values) / len(values)
        sorted_v = sorted(values)
        median = sorted_v[len(sorted_v) // 2]
        std = (sum((v - mean) ** 2 for v in values) / len(values)) ** 0.5
        print(f"{name:<15} {mean:>8.4f} {median:>8.4f} {std:>8.4f} {min(values):>8.4f} {max(values):>8.4f}")
    
    # Distribution of F1 scores
    print(f"\n--- F1 Score Distribution ---")
    bins = [(0.0, 0.0, 'F1 = 0 (complete miss)'),
            (0.001, 0.25, '0 < F1 ≤ 0.25'),
            (0.25, 0.5, '0.25 < F1 ≤ 0.5'),
            (0.5, 0.75, '0.5 < F1 ≤ 0.75'),
            (0.75, 1.0, '0.75 < F1 ≤ 1.0'),
            (1.0, 1.0, 'F1 = 1.0 (perfect)')]
    
    for lo, hi, label in bins:
        if lo == hi == 0.0:
            count = sum(1 for v in f1s if v == 0.0)
        elif lo == hi == 1.0:
            count = sum(1 for v in f1s if v == 1.0)
        else:
            count = sum(1 for v in f1s if lo < v <= hi)
        bar = '█' * count
        print(f"  {label:<30} {count:>4}  {bar}")
    
    # Perfect/zero breakdown
    perfect = sum(1 for f in f1s if f == 1.0)
    zero = sum(1 for f in f1s if f == 0.0)
    nonzero = sum(1 for f in f1s if f > 0.0)
    print(f"\n  Perfect predictions (F1=1.0): {perfect}/{n} ({perfect/n*100:.1f}%)")
    print(f"  Non-zero F1:                  {nonzero}/{n} ({nonzero/n*100:.1f}%)")
    print(f"  Complete misses (F1=0):        {zero}/{n} ({zero/n*100:.1f}%)")
    
    # Top performers and worst performers
    succeeded_with_f1 = [(e, float(e['f1'])) for e in succeeded 
                         if e.get('f1') and e['f1'] != '']
    succeeded_with_f1.sort(key=lambda x: x[1], reverse=True)
    
    print(f"\n--- Top 10 Best ---")
    print(f"  {'PDB':<8} {'MCSA':<8} {'F1':>6} {'Prec':>6} {'Rec':>6} {'Pred':>5} {'True':>5}")
    for e, f1 in succeeded_with_f1[:10]:
        print(f"  {e['pdb_id']:<8} {e['mcsa_id']:<8} {float(e['f1']):>6.3f} "
              f"{float(e['precision']):>6.3f} {float(e['recall']):>6.3f} "
              f"{e['n_predicted']:>5} {e['n_true']:>5}")
    
    print(f"\n--- Top 10 Worst ---")
    for e, f1 in succeeded_with_f1[-10:]:
        print(f"  {e['pdb_id']:<8} {e['mcsa_id']:<8} {float(e['f1']):>6.3f} "
              f"{float(e['precision']):>6.3f} {float(e['recall']):>6.3f} "
              f"{e['n_predicted']:>5} {e['n_true']:>5}")
    
    # Failure analysis
    if failed:
        print(f"\n--- Failures ---")
        error_counts = defaultdict(int)
        for e in failed:
            error = e.get('error', 'unknown').strip()[:60]
            error_counts[error] += 1
        for error, count in sorted(error_counts.items(), key=lambda x: -x[1]):
            print(f"  {count:>3}x  {error}")
    
    print("\n" + "=" * 70)
    
    # Save detailed JSON
    if args.output:
        results = {
            'n_total': total,
            'n_succeeded': len(succeeded),
            'n_failed': len(failed),
            'metrics': {
                'precision': {'mean': sum(precisions)/n, 'values': precisions},
                'recall': {'mean': sum(recalls)/n, 'values': recalls},
                'f1': {'mean': sum(f1s)/n, 'values': f1s},
            },
            'entries': [{'mcsa_id': e['mcsa_id'], 'pdb_id': e['pdb_id'],
                         'status': e['status'],
                         'precision': float(e['precision']) if e.get('precision') else None,
                         'recall': float(e['recall']) if e.get('recall') else None,
                         'f1': float(e['f1']) if e.get('f1') else None}
                        for e in entries]
        }
        with open(args.output, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"Saved detailed results to: {args.output}")


if __name__ == '__main__':
    main()
