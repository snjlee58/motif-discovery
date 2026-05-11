#!/usr/bin/env python3
"""
Compare F1 cohorts (perfect / high / low / zero) at the MSA-quality level to
identify whether the bimodality is driven by MSA depth/quality or by something
intrinsic to the catalytic columns in the failure cohort.

For each PDB:
  - Classify by F1 @ 1x:
      perfect: F1 == 1.0
      high   : 0.5 <= F1 < 1.0
      low    : 0    <  F1 < 0.5
      zero   : F1 == 0.0
  - MSA-level stats (from <pdb>_conservation.json statistics block):
      n_sequences, alignment_length, mean_conservation, mean_identity,
      mean_gap_frequency
  - Catalytic-column stats (mean over columns that map to M-CSA ground-truth):
      conservation, identity, gap_frequency, 3di_conservation
  - "Lift" = catalytic mean − MSA mean. Positive lift = catalytic columns
    are MORE conserved than typical for that MSA.

Interpretation matrix:

  - High lift in wins, zero lift in fails:
      conservation IS picking out catalytic residues, but only in some PDBs.
      Failures have MSAs where catalytic residues sit at high-divergence
      columns (functional but variable). 3Di redundancy is consistent with
      this — same columns are also variable in 3Di.
  - Similar lift in wins and fails:
      lift isn't the differentiator. Look at MSA depth (n_sequences) or
      other PDB-level factors.
  - Wins have much deeper MSAs than fails:
      cluster quality is the bottleneck — upstream issue.

Usage:
    python3 analysis/diagnose_msa_quality.py <rescore_dir> [--batch-dir <dir>]

Stdlib only.
"""

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path
from statistics import mean


COHORTS = ('perfect', 'high', 'low', 'zero')


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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('rescore_dir')
    ap.add_argument('--batch-dir', default=None)
    ap.add_argument('--out-tsv', default='analysis/msa_quality_table.tsv')
    args = ap.parse_args()

    rescore_dir = Path(args.rescore_dir)
    batch_dir = Path(args.batch_dir) if args.batch_dir else rescore_dir.parent
    if not rescore_dir.is_dir() or not batch_dir.is_dir():
        sys.exit(f"ERROR: directory missing  rescore={rescore_dir}  batch={batch_dir}")

    rows = []
    skipped = 0

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

        # Resid → columns
        resid_to_cols = defaultdict(set)
        for col, resid in mapping.items():
            resid_to_cols[resid].add(col)

        # Collect catalytic-column stats (one value per ground-truth residue,
        # using the best-conservation column if a residue maps to multiple cols).
        cat_cons, cat_ident, cat_gap, cat_3di = [], [], [], []
        for resid in gt:
            best = None
            for col in resid_to_cols.get(resid, ()):
                pos = resolve_pos(col, positions)
                if pos is None:
                    continue
                pd = positions[pos]
                c = pd.get('conservation') or 0.0
                if best is None or c > (best.get('conservation') or 0.0):
                    best = pd
            if best is None:
                continue
            cat_cons.append(best.get('conservation') or 0.0)
            cat_ident.append(best.get('identity') or 0.0)
            cat_gap.append(best.get('gap_frequency') or 0.0)
            t = best.get('3di_conservation')
            if t is not None:
                cat_3di.append(t)

        def safe_mean(xs):
            return mean(xs) if xs else None

        cm = safe_mean(cat_cons)
        ci = safe_mean(cat_ident)
        cg = safe_mean(cat_gap)
        c3 = safe_mean(cat_3di)

        row = {
            'pdb_id': pdb_id,
            'f1_at_1x': f1,
            'cohort': cohort,
            'n_true': n_true,
            'msa_n_sequences': cons.get('n_sequences'),
            'msa_alignment_length': cons.get('alignment_length'),
            'msa_mean_conservation': stats.get('mean_conservation'),
            'msa_mean_identity': stats.get('mean_identity'),
            'msa_mean_gap_frequency': stats.get('mean_gap_frequency'),
            'cat_mean_conservation': cm,
            'cat_mean_identity': ci,
            'cat_mean_gap_frequency': cg,
            'cat_mean_3di_conservation': c3,
            'lift_conservation': (cm - stats['mean_conservation']) if cm is not None and stats.get('mean_conservation') is not None else None,
            'lift_identity': (ci - stats['mean_identity']) if ci is not None and stats.get('mean_identity') is not None else None,
            'lift_gap_frequency': (cg - stats['mean_gap_frequency']) if cg is not None and stats.get('mean_gap_frequency') is not None else None,
        }
        rows.append(row)

    print(f"Loaded {len(rows)} PDBs from {rescore_dir}  (skipped {skipped})")
    cohort_counts = {c: sum(1 for r in rows if r['cohort'] == c) for c in COHORTS}
    print(f"Cohort sizes: " + ", ".join(f"{c}={cohort_counts[c]}" for c in COHORTS))
    print()

    def agg(cohort, key, decimals=3):
        vals = [r[key] for r in rows if r['cohort'] == cohort and r.get(key) is not None]
        if not vals:
            return None
        m = mean(vals)
        return round(m, decimals)

    def fmt(v, decimals):
        if v is None:
            return '—'
        if decimals == 0:
            return f'{v:.0f}'
        return f'{v:+.{decimals}f}' if decimals == 3 and isinstance(v, float) and 'lift' in '' else f'{v:.{decimals}f}'

    def print_section(title, fields):
        print(title)
        print("=" * 76)
        header = f"  {'stat':<32}" + "".join(f"{c:>10}" for c in COHORTS)
        print(header)
        for key, decimals, signed in fields:
            line = f"  {key:<32}"
            for c in COHORTS:
                v = agg(c, key, decimals)
                if v is None:
                    s = '—'
                else:
                    if signed:
                        s = f"{v:+.{decimals}f}"
                    elif decimals == 0:
                        s = f"{v:.0f}"
                    else:
                        s = f"{v:.{decimals}f}"
                line += f"{s:>10}"
            print(line)
        print()

    print_section("MSA-level stats (mean across PDBs in cohort)", [
        ('msa_n_sequences', 0, False),
        ('msa_alignment_length', 0, False),
        ('msa_mean_conservation', 3, False),
        ('msa_mean_identity', 3, False),
        ('msa_mean_gap_frequency', 3, False),
    ])
    print_section("Catalytic-column stats (mean across catalytic residues in cohort)", [
        ('cat_mean_conservation', 3, False),
        ('cat_mean_identity', 3, False),
        ('cat_mean_gap_frequency', 3, False),
        ('cat_mean_3di_conservation', 3, False),
    ])
    print_section("Lift (catalytic − MSA baseline; > 0 = catalytic is MORE conserved)", [
        ('lift_conservation', 3, True),
        ('lift_identity', 3, True),
        ('lift_gap_frequency', 3, True),
    ])

    # Per-PDB TSV
    out_path = Path(args.out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if rows:
        with open(out_path, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter='\t')
            w.writeheader()
            for r in rows:
                w.writerow({k: ('' if v is None else v) for k, v in r.items()})
        print(f"Wrote {out_path}")

    # Recommendation
    perfect_lift = agg('perfect', 'lift_conservation')
    zero_lift = agg('zero', 'lift_conservation')
    perfect_nseq = agg('perfect', 'msa_n_sequences', 0)
    zero_nseq = agg('zero', 'msa_n_sequences', 0)
    perfect_gap = agg('perfect', 'msa_mean_gap_frequency')
    zero_gap = agg('zero', 'msa_mean_gap_frequency')

    print()
    print("=" * 76)
    print("INTERPRETATION")
    print("=" * 76)
    if perfect_lift is not None and zero_lift is not None:
        dl = perfect_lift - zero_lift
        if dl > 0.05:
            print(f"  Catalytic-column conservation lift: perfect={perfect_lift:+.3f} vs zero={zero_lift:+.3f}")
            print(f"  → Wins have MUCH higher catalytic-column conservation than their MSA baseline,")
            print(f"    fails don't. The conservation signal works WHEN catalytic residues happen")
            print(f"    to be sequence-conserved. Failures are PDBs where catalytic residues sit at")
            print(f"    variable columns. 3Di redundancy explained: same columns are also variable")
            print(f"    in 3Di. Pivot: a different signal is needed for these (binding-pocket")
            print(f"    geometry, ML-learned scoring, or M-CSA-supervised weights).")
        elif abs(dl) < 0.02:
            print(f"  Catalytic-column conservation lift: perfect={perfect_lift:+.3f}, zero={zero_lift:+.3f}")
            print(f"  → Lifts are similar. The conservation signal differentiates catalytic residues")
            print(f"    equally well in wins and fails. The difference must be elsewhere — check")
            print(f"    MSA depth and gap stats above.")
        else:
            print(f"  Catalytic-column conservation lift: perfect={perfect_lift:+.3f}, zero={zero_lift:+.3f}")
            print(f"  → Modest difference (Δ={dl:+.3f}). Mixed signal; inspect table above by hand.")

    if perfect_nseq and zero_nseq:
        ratio = perfect_nseq / zero_nseq
        print()
        if ratio > 1.5:
            print(f"  MSA depth: perfect avg n_sequences = {perfect_nseq:.0f} vs zero = {zero_nseq:.0f} "
                  f"(ratio {ratio:.2f}×)")
            print(f"  → Wins have deeper MSAs. Cluster quality / AFDB membership is a factor.")
        elif ratio < 0.7:
            print(f"  MSA depth: perfect avg n_sequences = {perfect_nseq:.0f} vs zero = {zero_nseq:.0f}")
            print(f"  → Fails have MORE sequences than wins — counterintuitive, depth isn't the issue.")
            print(f"    Could be that deeper MSAs are NOISIER for catalytic residue prediction here.")
        else:
            print(f"  MSA depth: perfect avg n_sequences = {perfect_nseq:.0f} vs zero = {zero_nseq:.0f} "
                  f"(similar; depth not differentiating).")

    if perfect_gap is not None and zero_gap is not None:
        dg = zero_gap - perfect_gap
        if abs(dg) > 0.05:
            sign = "MORE" if dg > 0 else "FEWER"
            print(f"  MSA gappiness: zero cohort has {sign} gaps (Δ mean_gap_frequency = {dg:+.3f}).")

    print("=" * 76)


if __name__ == '__main__':
    main()
