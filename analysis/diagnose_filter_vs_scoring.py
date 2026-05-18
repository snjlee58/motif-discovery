#!/usr/bin/env python3
"""
For the F1=0-at-3x cohort whose catalytic residues ARE mapped to PDB residues
(the 'scoring_fail' cohort from diagnose_missing_residues.py), classify each
catalytic column by which step of the pipeline drops it:

  1. filter_gaps          — gap_frequency > 0.5 → --exclude-gaps drops it
  2. filter_identity      — identity < 0.2     → --min-identity 0.2 drops it
  3. pass_low_aa_low_3di    — passes filters; the column is NOT in top 3×n_true
                              by AA conservation alone, AND NOT in top 3×n_true
                              by 3Di conservation alone. Neither signal picks it up.
  4. pass_low_aa_high_3di   — passes filters; NOT in top 3×n_true by AA conservation,
                              BUT IS in top 3×n_true by 3Di conservation. These are
                              the residues that the user's 3Di-thesis would rescue
                              if 3Di were given more weight in the scoring formula.
  5. pass_high_conservation — passes filters, AND in top 3×n_true by AA conservation
                              alone, but the full scoring formula
                              (conservation + P2Rank + propensity + 3Di + clustering)
                              still didn't rank it in the top 3×n_true predicted set.
                              That means the OTHER features are actively downranking
                              the catalytic column — a scoring-weights issue.

Usage:
    python3 analysis/diagnose_filter_vs_scoring.py <rescore_dir> [--batch-dir <dir>]
        [--gap-threshold 0.5] [--identity-threshold 0.2]

Outputs:
    stdout — column-level cohort counts + sample PDBs per cohort
    analysis/filter_vs_scoring_table.tsv — one row per catalytic column in the cohort

Stdlib only.
"""

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path


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
    if 'metrics_by_top_n' in d and d['metrics_by_top_n']:
        labels = list(d['metrics_by_top_n'].keys())
        largest = labels[-1]
        pred_raw = d['metrics_by_top_n'][largest].get('predicted', []) or []
        top_n = d['metrics_by_top_n'][largest].get('top_n')
    else:
        largest = '1x'
        pred_raw = d.get('predicted', []) or []
        top_n = d.get('top_n')
    pred = {normalize_resid(v) for v in pred_raw}
    pred.discard(None)
    n_true = d.get('n_true') or (len(gt) if gt else None)
    return pdb_id, gt, pred, largest, top_n, n_true


def load_mapping(path):
    with open(path) as f:
        d = json.load(f)
    sub = d.get('mapping') if isinstance(d, dict) else None
    if not isinstance(sub, dict):
        return None
    # Returns {column_int: resid_normalized}
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
        d = json.load(f)
    positions = d.get('positions') or []
    # Return {position_1based: {gap_frequency, identity, conservation, consensus, ...}}
    out = {}
    for p in positions:
        pos = p.get('position')
        if pos is None:
            continue
        out[pos] = p
    return out


def find_pdb_dir(batch_dir, pdb_id):
    for variant in (pdb_id, pdb_id.upper(), pdb_id.lower()):
        cand = batch_dir / variant
        if cand.is_dir():
            return cand
    return None


def find_conservation_json(pdb_dir):
    # File pattern: <pdb_lower>_conservation.json
    candidates = list(pdb_dir.glob('*_conservation.json'))
    return candidates[0] if candidates else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('rescore_dir')
    ap.add_argument('--batch-dir', default=None)
    ap.add_argument('--gap-threshold', type=float, default=0.5,
                    help='--exclude-gaps drops columns with gap_frequency > this')
    ap.add_argument('--identity-threshold', type=float, default=0.2,
                    help='--min-identity drops columns with identity < this')
    ap.add_argument('--out-tsv', default='analysis/filter_vs_scoring_table.tsv')
    args = ap.parse_args()

    rescore_dir = Path(args.rescore_dir)
    if not rescore_dir.is_dir():
        sys.exit(f"ERROR: not a directory: {rescore_dir}")
    batch_dir = Path(args.batch_dir) if args.batch_dir else rescore_dir.parent
    if not batch_dir.is_dir():
        sys.exit(f"ERROR: batch dir not found: {batch_dir}")

    cohort_counts = Counter()        # per-residue cohort
    pdb_cohort_examples = {}         # cohort → list of (pdb_id, n_in_cohort, n_gt)
    table_rows = []
    skipped = Counter()
    largest_seen = None
    n_pdbs_in_cohort = 0

    for jf in sorted(rescore_dir.glob('*.json')):
        if jf.name == 'benchmark_summary.json':
            continue
        try:
            pdb_id, gt, pred, largest, top_n, n_true = load_rescore(jf)
        except (json.JSONDecodeError, KeyError):
            continue
        largest_seen = largest
        if not gt:
            continue
        # Only look at F1=0-at-largest-N (no overlap)
        if gt & pred:
            continue
        # Locate the batch sub-dir and required JSONs
        pdb_dir = find_pdb_dir(batch_dir, pdb_id)
        if pdb_dir is None:
            skipped['no_pdb_dir'] += 1
            continue
        mapping_path = pdb_dir / 'alignment_mapping.json'
        cons_path = find_conservation_json(pdb_dir)
        if not mapping_path.exists() or cons_path is None:
            skipped['missing_intermediates'] += 1
            continue
        try:
            mapping = load_mapping(mapping_path)  # {col_int: resid}
            positions = load_conservation(cons_path)  # {position_1based: dict}
        except json.JSONDecodeError:
            skipped['json_decode'] += 1
            continue
        if mapping is None:
            skipped['unknown_mapping_schema'] += 1
            continue

        # Reverse: resid → set of columns (a residue could appear at multiple cols
        # if the alignment has gaps that come back; treat each occurrence as one).
        resid_to_cols = {}
        for col, resid in mapping.items():
            resid_to_cols.setdefault(resid, set()).add(col)

        # Only the 'scoring_fail' cohort matters here: residues that ARE mapped.
        gt_mapped = [r for r in gt if r in resid_to_cols]
        if not gt_mapped:
            # mapping_fail_all — handled by other diagnostic
            continue

        # Rank columns by AA conservation and by 3Di conservation independently.
        # Columns missing a 3Di score are treated as score 0 for ranking purposes.
        cons_sorted = sorted(
            ((pos, pdata.get('conservation', 0.0) or 0.0) for pos, pdata in positions.items()),
            key=lambda x: -x[1],
        )
        cons_rank = {pos: i + 1 for i, (pos, _) in enumerate(cons_sorted)}

        cons3di_sorted = sorted(
            ((pos, pdata.get('3di_conservation', 0.0) or 0.0) for pos, pdata in positions.items()),
            key=lambda x: -x[1],
        )
        cons3di_rank = {pos: i + 1 for i, (pos, _) in enumerate(cons3di_sorted)}
        # Detect whether 3Di was actually computed (else all zeros → ranks are arbitrary)
        any_3di = any((pdata.get('3di_conservation') or 0) > 0 for pdata in positions.values())

        # how many "top" by conservation alone we'd accept
        top_k = top_n or (len(gt) * 3)  # fallback to 3x if not set

        # Mapping uses column indices that are 1-based to match position numbering
        # (map_alignment_to_pdb.py preserves the alignment column index;
        # score_conservation.py uses position = col_idx + 1). We'll try the
        # column-as-is first, then column+1, and use whichever matches a position.
        def resolve_position(col):
            if col in positions:
                return col
            if (col + 1) in positions:
                return col + 1
            return None

        n_in_cohort = len(gt_mapped)
        n_pdbs_in_cohort += 1
        per_pdb_counter = Counter()
        for resid in gt_mapped:
            cols = resid_to_cols[resid]
            # Use the best (lowest cohort) column for this residue.
            best_cohort = None
            best_row = None
            for col in cols:
                pos = resolve_position(col)
                if pos is None:
                    cohort = 'position_lookup_failed'
                    row = {
                        'pdb_id': pdb_id, 'resid': resid, 'column': col, 'position': '',
                        'gap_frequency': '', 'identity': '', 'conservation': '',
                        'conservation_rank': '', 'top_k': top_k, 'cohort': cohort,
                    }
                else:
                    pdata = positions[pos]
                    gf = pdata.get('gap_frequency') or 0.0
                    ident = pdata.get('identity') or 0.0
                    cons = pdata.get('conservation') or 0.0
                    cons3di = pdata.get('3di_conservation')
                    cons3di_val = cons3di if cons3di is not None else 0.0
                    rank = cons_rank.get(pos, -1)
                    rank_3di = cons3di_rank.get(pos, -1) if any_3di else -1
                    if gf > args.gap_threshold:
                        cohort = 'filter_gaps'
                    elif ident < args.identity_threshold:
                        cohort = 'filter_identity'
                    elif rank <= top_k:
                        cohort = 'pass_high_conservation'
                    elif any_3di and 0 < rank_3di <= top_k:
                        cohort = 'pass_low_aa_high_3di'
                    else:
                        cohort = 'pass_low_aa_low_3di'
                    row = {
                        'pdb_id': pdb_id, 'resid': resid, 'column': col, 'position': pos,
                        'gap_frequency': f'{gf:.3f}', 'identity': f'{ident:.3f}',
                        'conservation': f'{cons:.3f}',
                        'conservation_rank': rank,
                        '3di_conservation': f'{cons3di_val:.3f}' if cons3di is not None else '',
                        '3di_rank': rank_3di if any_3di else '',
                        'top_k': top_k, 'cohort': cohort,
                    }
                # Track the best (lowest-severity) cohort if multiple columns map to same resid
                severity = {'filter_gaps': 4, 'filter_identity': 3,
                            'pass_low_aa_low_3di': 2, 'pass_low_aa_high_3di': 1,
                            'pass_high_conservation': 0,
                            'position_lookup_failed': 5}
                if best_cohort is None or severity.get(cohort, 5) < severity.get(best_cohort, 5):
                    best_cohort = cohort
                    best_row = row
            cohort_counts[best_cohort] += 1
            per_pdb_counter[best_cohort] += 1
            table_rows.append(best_row)

        # For per-PDB example tracking, label by the dominant cohort
        dominant = per_pdb_counter.most_common(1)[0][0]
        pdb_cohort_examples.setdefault(dominant, []).append((pdb_id, per_pdb_counter[dominant], n_in_cohort))

    total_resids = sum(cohort_counts.values())
    print(f"Rescore dir: {rescore_dir}")
    print(f"Batch dir:   {batch_dir}")
    print(f"Largest top-N label: {largest_seen}")
    print(f"Filters being modeled: gap_freq > {args.gap_threshold}, identity < {args.identity_threshold}")
    print()
    print(f"Catalytic residues analyzed (F1=0-at-largest-N AND mapped): {total_resids} across {n_pdbs_in_cohort} PDBs")
    if skipped:
        print(f"Skipped PDBs: {dict(skipped)}")
    print()
    print("=" * 70)
    print(f"  {'cohort':<28}{'count':>8}{'%':>8}")
    print("=" * 70)
    order = ['filter_gaps', 'filter_identity',
             'pass_low_aa_low_3di', 'pass_low_aa_high_3di',
             'pass_high_conservation', 'position_lookup_failed']
    for c in order:
        n = cohort_counts.get(c, 0)
        pct = 100 * n / total_resids if total_resids else 0
        print(f"  {c:<28}{n:>8}{pct:>7.1f}%")

    print()
    print("Per-PDB dominant cohort (PDBs grouped by their majority residue cohort):")
    for c in order:
        pdbs = pdb_cohort_examples.get(c, [])
        if pdbs:
            sample = [p for p, _, _ in pdbs[:15]]
            print(f"  {c} ({len(pdbs)} PDBs): {sample}{' …' if len(pdbs) > 15 else ''}")

    out_path = Path(args.out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=[
            'pdb_id', 'resid', 'column', 'position', 'gap_frequency', 'identity',
            'conservation', 'conservation_rank',
            '3di_conservation', '3di_rank',
            'top_k', 'cohort'], delimiter='\t')
        w.writeheader()
        w.writerows(table_rows)
    print(f"\nWrote {out_path}  ({len(table_rows)} rows)")

    print()
    print("=" * 70)
    if total_resids == 0:
        print("No residues to analyze. Check rescore_dir / batch_dir paths.")
    else:
        n_filter = cohort_counts.get('filter_gaps', 0) + cohort_counts.get('filter_identity', 0)
        n_3di_rescue = cohort_counts.get('pass_low_aa_high_3di', 0)
        n_pass_low = cohort_counts.get('pass_low_aa_low_3di', 0)
        n_pass_high = cohort_counts.get('pass_high_conservation', 0)

        print(f"RECOMMENDATION SUMMARY ({total_resids} catalytic residues in F1=0 cohort):")
        print()
        if n_filter:
            print(f"  • {n_filter} ({100*n_filter/total_resids:.0f}%) lost to filters → relax "
                  f"--exclude-gaps and/or --min-identity (cheapest fix).")
        if n_3di_rescue:
            print(f"  • {n_3di_rescue} ({100*n_3di_rescue/total_resids:.0f}%) rescuable by 3Di "
                  f"conservation alone → up-weight 3Di in scoring (validates the thesis).")
        if n_pass_low:
            print(f"  • {n_pass_low} ({100*n_pass_low/total_resids:.0f}%) miss in BOTH AA & 3Di → "
                  f"neither signal works for these. Hardest cohort. Likely needs another signal "
                  f"(P2Rank, propensity, cluster context) or these are M-CSA-label edge cases.")
        if n_pass_high:
            print(f"  • {n_pass_high} ({100*n_pass_high/total_resids:.0f}%) reached top-N by "
                  f"AA conservation but were downranked by other features → ablate or re-tune weights.")

        def has_meaningful_3di(rows):
            for r in rows:
                s = r.get('3di_conservation')
                if not s:
                    continue
                try:
                    if float(s) > 0:
                        return True
                except (TypeError, ValueError):
                    continue
            return False

        if not has_meaningful_3di(table_rows):
            print()
            print("  WARNING: 3di_conservation is empty/zero throughout. Either the pipeline")
            print("  didn't save 3Di scores for this batch, or 3di_conservation was computed")
            print("  but all values are 0. Re-run pipeline with the 3Di alignment passed to")
            print("  score_conservation.py to populate this signal — then re-run this diagnostic.")
    print("=" * 70)


if __name__ == '__main__':
    main()
