#!/usr/bin/env python3
"""
Classify the F1=0-at-3x failure cohort by where in the pipeline the catalytic
residue is being lost.

For each PDB whose rescore JSON shows F1=0 even at the largest top-N multiplier,
this checks whether the M-CSA ground-truth residue IDs are present in that PDB's
`alignment_mapping.json`:

  - NOT present in mapping        → map_alignment_to_pdb.py dropped the residue
                                    (M-CSA numbering vs. PDB auth_resid mismatch,
                                    chain selection, insertion codes, etc.)
  - Present in mapping but missed → the column WAS mapped, so the residue ID is
                                    available; conservation scoring + filters
                                    didn't rank it in the top 3N. Either the
                                    --exclude-gaps / --min-identity 0.2 filter
                                    excluded the column, or its conservation
                                    score was too low.

Usage:
    python3 analysis/diagnose_missing_residues.py <rescore_dir> [--batch-dir <dir>]

If --batch-dir is omitted, the script assumes the rescore dir is a subdirectory
of the original batch dir (the standard layout from rescore_batch.sh).

Outputs:
    stdout — per-cohort counts + sample PDB lists
    analysis/missing_residues_table.tsv — one row per F1=0-at-largest-N PDB

Stdlib only.
"""

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path


def load_rescore(jf):
    """Return (pdb_id, ground_truth_set, predicted_at_largest_N, largest_label)."""
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
    else:
        largest = '1x'
        pred_raw = d.get('predicted', []) or []

    pred = {normalize_resid(v) for v in pred_raw}
    pred.discard(None)
    return pdb_id, gt, pred, largest


def normalize_resid(v):
    """Residue IDs may be int or numeric string. Insertion-code strings (e.g. "70A")
    are kept as strings."""
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
        return s  # insertion-coded resid like "70A"


def load_mapping(path):
    """Return set of all residue IDs that were mapped in alignment_mapping.json.

    Schema varies — handle common shapes defensively:
      {col: resid, ...}
      {"mapping": {col: resid, ...}, ...metadata}
      {"alignment_to_pdb": {col: resid, ...}}
    """
    with open(path) as f:
        d = json.load(f)
    candidate_dicts = []
    if isinstance(d, dict):
        # Direct mapping at top level (keys are columns, values are resids)
        if all(isinstance(v, (int, str, float, type(None))) for v in d.values()):
            candidate_dicts.append(d)
        for key in ('mapping', 'alignment_to_pdb', 'column_to_resid', 'col_to_resid'):
            sub = d.get(key)
            if isinstance(sub, dict):
                candidate_dicts.append(sub)
    if not candidate_dicts:
        return None  # unrecognized schema
    mapped = set()
    for sub in candidate_dicts:
        for v in sub.values():
            r = normalize_resid(v)
            if r is not None:
                mapped.add(r)
    return mapped


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('rescore_dir')
    ap.add_argument('--batch-dir', default=None,
                    help='Original batch dir containing <PDB>/alignment_mapping.json. '
                         'Defaults to rescore_dir.parent.')
    ap.add_argument('--out-tsv', default='analysis/missing_residues_table.tsv')
    args = ap.parse_args()

    rescore_dir = Path(args.rescore_dir)
    if not rescore_dir.is_dir():
        sys.exit(f"ERROR: not a directory: {rescore_dir}")
    batch_dir = Path(args.batch_dir) if args.batch_dir else rescore_dir.parent
    if not batch_dir.is_dir():
        sys.exit(f"ERROR: batch dir not found: {batch_dir}")

    cohorts = Counter()
    examples = {
        'mapping_fail_all': [],   # 0/N ground-truth residues in alignment_mapping
        'mapping_fail_partial': [],  # some missing, some mapped
        'scoring_fail': [],       # all N mapped, but none predicted at largest N
        'no_mapping_file': [],
        'mapping_schema_unknown': [],
        'no_ground_truth': [],
    }
    table_rows = []
    largest_label_seen = None

    for jf in sorted(rescore_dir.glob('*.json')):
        if jf.name == 'benchmark_summary.json':
            continue
        try:
            pdb_id, gt, pred, largest = load_rescore(jf)
        except (json.JSONDecodeError, KeyError):
            continue
        largest_label_seen = largest

        if not gt:
            examples['no_ground_truth'].append(pdb_id)
            continue

        # Only interested in F1=0 at largest N
        if gt & pred:
            continue

        # Find the alignment_mapping.json for this PDB. Try common case variants.
        mapping_path = None
        for variant in (pdb_id, pdb_id.upper(), pdb_id.lower()):
            cand = batch_dir / variant / 'alignment_mapping.json'
            if cand.exists():
                mapping_path = cand
                break
        if mapping_path is None:
            cohorts['no_mapping_file'] += 1
            examples['no_mapping_file'].append(pdb_id)
            table_rows.append({
                'pdb_id': pdb_id, 'cohort': 'no_mapping_file',
                'n_ground_truth': len(gt), 'n_mapped': 0, 'n_in_predicted_at_largest_N': 0,
                'largest_top_n_label': largest,
                'ground_truth': ','.join(str(x) for x in sorted(gt, key=str)),
                'note': 'alignment_mapping.json not found under batch_dir',
            })
            continue

        mapped = load_mapping(mapping_path)
        if mapped is None:
            cohorts['mapping_schema_unknown'] += 1
            examples['mapping_schema_unknown'].append(pdb_id)
            table_rows.append({
                'pdb_id': pdb_id, 'cohort': 'mapping_schema_unknown',
                'n_ground_truth': len(gt), 'n_mapped': 0, 'n_in_predicted_at_largest_N': 0,
                'largest_top_n_label': largest,
                'ground_truth': ','.join(str(x) for x in sorted(gt, key=str)),
                'note': f'unrecognized schema in {mapping_path}',
            })
            continue

        n_gt = len(gt)
        n_mapped_gt = sum(1 for r in gt if r in mapped)

        if n_mapped_gt == 0:
            cohort = 'mapping_fail_all'
        elif n_mapped_gt < n_gt:
            cohort = 'mapping_fail_partial'
        else:
            cohort = 'scoring_fail'

        cohorts[cohort] += 1
        examples[cohort].append(pdb_id)
        table_rows.append({
            'pdb_id': pdb_id, 'cohort': cohort,
            'n_ground_truth': n_gt,
            'n_mapped': n_mapped_gt,
            'n_in_predicted_at_largest_N': 0,  # by construction (we filtered F1=0)
            'largest_top_n_label': largest,
            'ground_truth': ','.join(str(x) for x in sorted(gt, key=str)),
            'note': '',
        })

    total_f1_zero = sum(cohorts[k] for k in
                       ('mapping_fail_all', 'mapping_fail_partial', 'scoring_fail',
                        'no_mapping_file', 'mapping_schema_unknown'))

    print(f"Rescore dir: {rescore_dir}")
    print(f"Batch dir:   {batch_dir}")
    print(f"Largest top-N label compared: {largest_label_seen}")
    print()
    print(f"F1=0-at-{largest_label_seen} cohort: {total_f1_zero} PDBs")
    print("=" * 70)
    print(f"  mapping_fail_all     : {cohorts['mapping_fail_all']:>4}  "
          f"({100*cohorts['mapping_fail_all']/max(total_f1_zero,1):.1f}%)  "
          f"— 0/N ground-truth residues appear in alignment_mapping.json")
    print(f"  mapping_fail_partial : {cohorts['mapping_fail_partial']:>4}  "
          f"({100*cohorts['mapping_fail_partial']/max(total_f1_zero,1):.1f}%)  "
          f"— some missing, some mapped (mixed mapping issue)")
    print(f"  scoring_fail         : {cohorts['scoring_fail']:>4}  "
          f"({100*cohorts['scoring_fail']/max(total_f1_zero,1):.1f}%)  "
          f"— N/N mapped, but none predicted at top {largest_label_seen}")
    if cohorts['no_mapping_file']:
        print(f"  no_mapping_file      : {cohorts['no_mapping_file']:>4}  "
              f"— alignment_mapping.json missing from batch dir")
    if cohorts['mapping_schema_unknown']:
        print(f"  mapping_schema_unknown: {cohorts['mapping_schema_unknown']:>4}  "
              f"— mapping.json present but schema not recognized")

    def show(label, key, n=15):
        items = examples[key]
        if items:
            print(f"\n  {label} ({len(items)}): "
                  f"{items[:n]}{' …' if len(items) > n else ''}")

    show('mapping_fail_all',     'mapping_fail_all')
    show('mapping_fail_partial', 'mapping_fail_partial')
    show('scoring_fail',         'scoring_fail')
    show('no_mapping_file',      'no_mapping_file')
    show('mapping_schema_unknown','mapping_schema_unknown')
    if examples['no_ground_truth']:
        show('no_ground_truth (skipped)', 'no_ground_truth')

    # Write per-PDB table
    out_path = Path(args.out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ['pdb_id', 'cohort', 'n_ground_truth', 'n_mapped',
                  'n_in_predicted_at_largest_N', 'largest_top_n_label',
                  'ground_truth', 'note']
    with open(out_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        w.writeheader()
        w.writerows(table_rows)
    print(f"\nWrote {out_path}  ({len(table_rows)} rows — every F1=0-at-{largest_label_seen} PDB)")

    # Recommendation
    print()
    print("=" * 70)
    mf = cohorts['mapping_fail_all'] + cohorts['mapping_fail_partial']
    sf = cohorts['scoring_fail']
    if mf > sf * 1.5:
        print(f"RECOMMENDATION: mapping failures dominate ({mf} vs scoring {sf}).")
        print("  → Look at map_alignment_to_pdb.py — residue numbering / chain / insertion-code")
        print("    issues are likely losing catalytic residues before they can be scored.")
        print("  → Inspect a few mapping_fail_all PDBs by hand: open their alignment_mapping.json")
        print("    and compare the values set to the M-CSA ground-truth residue IDs.")
    elif sf > mf * 1.5:
        print(f"RECOMMENDATION: scoring failures dominate ({sf} vs mapping {mf}).")
        print("  → Residues ARE in alignment_mapping but conservation scoring isn't ranking them")
        print("    in the top 3N. Tune extract_top_conserved.py filters or scoring weights.")
        print("  → Check --exclude-gaps --min-identity 0.2 — are these excluding the columns")
        print("    that contain the catalytic residues?")
    else:
        print(f"RECOMMENDATION: split fairly even (mapping {mf} vs scoring {sf}).")
        print("  → Two parallel fixes needed. Tackle mapping first (lossy data) then scoring.")
    print("=" * 70)


if __name__ == '__main__':
    main()
