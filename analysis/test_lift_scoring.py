#!/usr/bin/env python3
"""
Re-rank conservation columns using lift-based scoring and compare new F1 to the
original conservation-based scoring across F1 cohorts.

Hypothesis (from diagnose_msa_quality.py): the F1=0 cohort has catalytic columns
with NEGATIVE lift_conservation (catalytic columns are LESS conserved than the
MSA baseline), even though lift_identity remains positive. Re-ranking by lift
over the MSA baseline should recover those residues.

Variants tested:
  lift_conservation   = column.conservation - msa.mean_conservation
  lift_identity       = column.identity     - msa.mean_identity
  lift_3di            = column.3di_conservation - mean(3di_conservation)
  lift_combined       = average of lift_identity + lift_3di (both signals lifted)

For each variant and each PDB:
  1. Read alignment_mapping.json → column → residue map
  2. Read <pdb>_conservation.json → per-column stats + MSA-level statistics
  3. Compute lift for every mapped column
  4. Take top n_true columns by that lift → predicted residue set
  5. Compute F1 vs M-CSA ground truth
  6. Aggregate by original-F1 cohort (perfect/high/low/zero) — does the zero
     cohort recover?

Output:
  stdout: per-cohort F1 (original vs each variant) + total
  analysis/lift_scoring_results.tsv: per-PDB original F1 + new F1 per variant

Stdlib only. Read-only — does NOT touch pipeline.sh or benchmark_mcsa.py.

Usage:
    python3 analysis/test_lift_scoring.py <rescore_dir> [--batch-dir <dir>]
"""

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path
from statistics import mean


COHORTS = ('perfect', 'high', 'low', 'zero')
VARIANTS = ('lift_conservation', 'lift_identity', 'lift_3di', 'lift_combined')


def normalize_resid(v):
    if v is None:
        return None
    if isinstance(v, int):
        return v
    s = str(v).strip()
    if s == '':
        return None
    try:
        return int(s)
    except ValueError:
        return s


def load_rescore(jf):
    with open(jf) as f:
        d = json.load(f)
    pdb_id = d.get('pdb_id')
    gt_raw = d.get('mcsa_ground_truth') or []
    gt = {normalize_resid(v) for v in gt_raw}
    gt.discard(None)
    if 'metrics_by_top_n' in d and '1x' in d['metrics_by_top_n']:
        m = d['metrics_by_top_n']['1x']['metrics']
    else:
        m = d.get('metrics', {})
    f1 = m.get('f1')
    n_true = d.get('n_true') or len(gt)
    return pdb_id, gt, f1, n_true


def load_mapping(path):
    with open(path) as f:
        d = json.load(f)
    sub = d.get('mapping') if isinstance(d, dict) else None
    if not isinstance(sub, dict):
        return None
    out = {}
    for col_str, resid in sub.items():
        try:
            col = int(col_str)
        except (ValueError, TypeError):
            continue
        r = normalize_resid(resid)
        if r is not None:
            out[col] = r
    return out


def load_conservation(path):
    with open(path) as f:
        return json.load(f)


def find_pdb_dir(batch_dir, pdb_id):
    for variant in (pdb_id, pdb_id.upper(), pdb_id.lower()):
        cand = batch_dir / variant
        if cand.is_dir():
            return cand
    return None


def find_conservation_json(pdb_dir):
    candidates = list(pdb_dir.glob('*_conservation.json'))
    return candidates[0] if candidates else None


def cohort_for(f1):
    if f1 is None:
        return None
    if f1 == 1.0:
        return 'perfect'
    if f1 >= 0.5:
        return 'high'
    if f1 > 0:
        return 'low'
    return 'zero'


def resolve_pos(col, positions):
    if col in positions:
        return col
    if (col + 1) in positions:
        return col + 1
    return None


def f1_score(predicted, ground_truth, n_predicted):
    """F1 where precision uses n_predicted (matches benchmark_mcsa.py convention)."""
    if not ground_truth or not n_predicted:
        return 0.0, 0, 0.0, 0.0
    tp = len(predicted & ground_truth)
    p = tp / n_predicted
    r = tp / len(ground_truth)
    f1 = (2 * p * r / (p + r)) if (p + r) else 0.0
    return f1, tp, p, r


def rerank_and_score(positions, mapping, gt, n_true, msa_mean_id, msa_mean_cons, score_fn):
    """For one PDB, rank mapped columns by score_fn and pick top n_true."""
    # Compute per-resid best score (a resid may map to multiple columns)
    resid_scores = {}
    for col, resid in mapping.items():
        pos = resolve_pos(col, positions)
        if pos is None:
            continue
        pdata = positions[pos]
        s = score_fn(pdata, msa_mean_id, msa_mean_cons)
        if s is None:
            continue
        # Take the best column score per residue
        if resid not in resid_scores or s > resid_scores[resid]:
            resid_scores[resid] = s
    if not resid_scores:
        return 0.0, 0, set()
    ranked = sorted(resid_scores.items(), key=lambda kv: -kv[1])
    predicted = {r for r, _ in ranked[:n_true]}
    f1, tp, _, _ = f1_score(predicted, gt, n_true)
    return f1, tp, predicted


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('rescore_dir')
    ap.add_argument('--batch-dir', default=None)
    ap.add_argument('--out-tsv', default='analysis/lift_scoring_results.tsv')
    args = ap.parse_args()

    rescore_dir = Path(args.rescore_dir)
    batch_dir = Path(args.batch_dir) if args.batch_dir else rescore_dir.parent
    if not rescore_dir.is_dir() or not batch_dir.is_dir():
        sys.exit(f"ERROR: directory missing  rescore={rescore_dir}  batch={batch_dir}")

    # Per-PDB rows
    rows = []
    skipped = 0

    # Score functions per variant.
    def sf_lift_cons(p, mid, mcons):
        c = p.get('conservation')
        return None if c is None or mcons is None else c - mcons

    def sf_lift_id(p, mid, mcons):
        i = p.get('identity')
        return None if i is None or mid is None else i - mid

    def sf_lift_3di(p, mid, mcons):
        t = p.get('3di_conservation')
        return None if t is None else t  # mean handled below by per-PDB centering

    def sf_lift_combined(p, mid, mcons):
        i = p.get('identity')
        t = p.get('3di_conservation')
        if i is None or mid is None or t is None:
            return None
        return (i - mid) + t  # simple sum; 3di is already inverted-entropy

    score_fns = {
        'lift_conservation': sf_lift_cons,
        'lift_identity':     sf_lift_id,
        'lift_3di':          sf_lift_3di,
        'lift_combined':     sf_lift_combined,
    }

    for jf in sorted(rescore_dir.glob('*.json')):
        if jf.name == 'benchmark_summary.json':
            continue
        try:
            pdb_id, gt, f1, n_true = load_rescore(jf)
        except (json.JSONDecodeError, KeyError):
            skipped += 1
            continue
        cohort = cohort_for(f1)
        if cohort is None or not gt:
            skipped += 1
            continue

        pdb_dir = find_pdb_dir(batch_dir, pdb_id)
        if pdb_dir is None:
            skipped += 1
            continue
        map_path = pdb_dir / 'alignment_mapping.json'
        cons_path = find_conservation_json(pdb_dir)
        if not map_path.exists() or cons_path is None:
            skipped += 1
            continue
        try:
            mapping = load_mapping(map_path)
            cons = load_conservation(cons_path)
        except json.JSONDecodeError:
            skipped += 1
            continue
        if mapping is None:
            skipped += 1
            continue

        positions = {p.get('position'): p for p in (cons.get('positions') or []) if p.get('position')}
        stats = cons.get('statistics') or {}
        msa_mean_id = stats.get('mean_identity')
        msa_mean_cons = stats.get('mean_conservation')

        # 3Di "lift" needs a per-MSA mean; compute it once here.
        all_3di = [p.get('3di_conservation') for p in positions.values()
                   if p.get('3di_conservation') is not None]
        msa_mean_3di = mean(all_3di) if all_3di else None

        # Patch the 3Di score fn for this PDB to subtract the local 3Di mean.
        def sf_lift_3di_local(p, mid, mcons, m3=msa_mean_3di):
            t = p.get('3di_conservation')
            return None if t is None or m3 is None else t - m3

        def sf_lift_combined_local(p, mid, mcons, m3=msa_mean_3di):
            i = p.get('identity')
            t = p.get('3di_conservation')
            if i is None or mid is None or t is None or m3 is None:
                return None
            return (i - mid) + (t - m3)

        score_fns_local = dict(score_fns)
        score_fns_local['lift_3di'] = sf_lift_3di_local
        score_fns_local['lift_combined'] = sf_lift_combined_local

        row = {
            'pdb_id': pdb_id,
            'cohort': cohort,
            'original_f1': f1,
            'n_true': n_true,
        }
        for vname, fn in score_fns_local.items():
            new_f1, tp, _ = rerank_and_score(positions, mapping, gt, n_true,
                                              msa_mean_id, msa_mean_cons, fn)
            row[f'{vname}_f1'] = round(new_f1, 4)
            row[f'{vname}_tp'] = tp
        rows.append(row)

    if not rows:
        sys.exit("No PDBs processed; check paths and JSON schemas.")

    # Aggregate
    print(f"Loaded {len(rows)} PDBs from {rescore_dir}  (skipped {skipped})")
    cohort_counts = {c: sum(1 for r in rows if r['cohort'] == c) for c in COHORTS}
    print(f"Cohort sizes: " + ", ".join(f"{c}={cohort_counts[c]}" for c in COHORTS))
    print()

    print("Mean F1 by cohort (original vs lift variants)")
    print("=" * 82)
    header = f"  {'cohort':<10}{'n':>5}{'orig':>10}"
    for v in VARIANTS:
        header += f"{v:>16}"
    print(header)
    for c in COHORTS:
        line = f"  {c:<10}{cohort_counts[c]:>5}"
        # Original
        f1s = [r['original_f1'] for r in rows if r['cohort'] == c and r['original_f1'] is not None]
        orig_mean = mean(f1s) if f1s else 0.0
        line += f"{orig_mean:>10.3f}"
        for v in VARIANTS:
            vs = [r[f'{v}_f1'] for r in rows if r['cohort'] == c]
            vm = mean(vs) if vs else 0.0
            delta = vm - orig_mean
            line += f"  {vm:>6.3f} ({delta:+.3f})"
        print(line)
    # Total
    line = f"  {'TOTAL':<10}{len(rows):>5}"
    f1s = [r['original_f1'] for r in rows if r['original_f1'] is not None]
    orig_mean = mean(f1s) if f1s else 0.0
    line += f"{orig_mean:>10.3f}"
    for v in VARIANTS:
        vs = [r[f'{v}_f1'] for r in rows]
        vm = mean(vs) if vs else 0.0
        delta = vm - orig_mean
        line += f"  {vm:>6.3f} ({delta:+.3f})"
    print(line)
    print()

    # F1=0 rescue table
    print("F1=0-cohort rescue counts (how many of the failure PDBs each variant lifts to F1>0)")
    print("=" * 82)
    zero_rows = [r for r in rows if r['cohort'] == 'zero']
    n_zero = len(zero_rows)
    if n_zero:
        print(f"  variant            rescued (F1>0)    rescued to F1≥0.5    perfect (F1=1)")
        for v in VARIANTS:
            rescued = sum(1 for r in zero_rows if r[f'{v}_f1'] > 0)
            strong = sum(1 for r in zero_rows if r[f'{v}_f1'] >= 0.5)
            perfect = sum(1 for r in zero_rows if r[f'{v}_f1'] >= 1.0)
            print(f"  {v:<18}{rescued:>10}/{n_zero}   {strong:>14}/{n_zero}     {perfect:>10}/{n_zero}")
    print()

    # Per-PDB TSV
    out_path = Path(args.out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with open(out_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {out_path}")

    # Recommendation
    print()
    print("=" * 82)
    print("INTERPRETATION")
    print("=" * 82)
    # Find best variant by total F1
    best_v = None
    best_delta = -1
    for v in VARIANTS:
        vs = [r[f'{v}_f1'] for r in rows]
        delta = mean(vs) - orig_mean
        if delta > best_delta:
            best_delta = delta
            best_v = v
    if best_delta > 0.05:
        print(f"  Best variant: {best_v}  (mean F1 lift {best_delta:+.3f} vs original).")
        print(f"  → This is the candidate scoring formula. Modify benchmark_mcsa.py to use")
        print(f"    {best_v}-style ranking and re-run rescore on a fresh batch to validate.")
    elif best_delta > 0.01:
        print(f"  Modest improvement: {best_v} mean F1 lift {best_delta:+.3f}. Worth pursuing")
        print(f"  but not a silver bullet — combine with other signals (P2Rank, propensity).")
    else:
        print(f"  No variant clearly beats original (best Δ={best_delta:+.3f} via {best_v}).")
        print(f"  Lift-based scoring alone is not enough — the deep-MSA failures may need")
        print(f"  a fundamentally different signal (binding-pocket geometry, ML scoring).")
    print("=" * 82)


if __name__ == '__main__':
    main()
