#!/usr/bin/env python3
"""
Diagnose bimodal F1 performance on the M-CSA benchmark — stdlib only.

Joins a `benchmark_summary.tsv` (from `benchmark/rescore_batch.sh`) with
`analysis/cluster_sizes.tsv` and looks for protein-/cluster-level properties
that differentiate the success cohort from the failure cohort.

Usage:
    python3 analysis/diagnose_bimodal.py <benchmark_summary.tsv> \\
        [--cluster-tsv analysis/cluster_sizes.tsv] \\
        [--out-dir analysis]

Outputs:
    analysis/diagnostic_table.tsv         joined per-protein wide table
    analysis/diagnostic_long.tsv          long-format (pdb_id, dimension, value, f1)
                                          for easy plotting in any tool
    stdout                                Spearman ρ, Mann-Whitney U, bucket means,
                                          ASCII F1 histogram, recommendation

Degrades gracefully — if cluster_sizes.tsv has only the basic columns
(mcsa_id, pdb_id, uniprot, rep_id, cluster_size), it stratifies on those.
Re-run analyze_cluster_sizes.py with AFDB file 2 to unlock richer dimensions.
"""

import argparse
import csv
import math
import sys
from pathlib import Path


NUMERIC_CANDIDATES = [
    'n_true',
    'cluster_size_afdb50', 'cluster_size_full', 'cluster_size',
    'dedup_ratio', 'avg_pLDDT', 'avg_seq_len', 'nMem_file2',
]
CATEGORICAL_CANDIDATES = ['isDark']


# ---------- IO ----------

def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(path, rows, fieldnames):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in fieldnames})


def as_float(v):
    """Coerce string → float, return None for blanks/junk."""
    if v is None:
        return None
    s = str(v).strip()
    if s == '' or s.upper() == 'NA':
        return None
    try:
        return float(s)
    except ValueError:
        return None


def cluster_lookup_ok(row):
    """A row has a real cluster lookup only if uniprot resolved and cluster_size > 0.

    cluster_sizes.tsv from analyze_cluster_sizes.py marks failed UniProt
    bulk-API resolutions as uniprot='NOT_FOUND' rep_id='NOT_FOUND' cluster_size=0.
    Those rows should not be used for cluster-size correlations.
    """
    uniprot = (row.get('uniprot') or '').strip()
    rep_id = (row.get('rep_id') or '').strip()
    cs = as_float(row.get('cluster_size'))
    cs_alt = as_float(row.get('cluster_size_afdb50')) if 'cluster_size_afdb50' in row else None
    cs_full = as_float(row.get('cluster_size_full')) if 'cluster_size_full' in row else None
    if uniprot == 'NOT_FOUND' or rep_id == 'NOT_FOUND':
        return False
    return any(v is not None and v > 0 for v in (cs, cs_alt, cs_full))


# ---------- Join ----------

def join(summary_rows, cluster_rows):
    """Left-join summary with cluster on uppercase pdb_id."""
    by_key = {r['pdb_id'].upper(): r for r in cluster_rows}
    cluster_cols = [c for c in (cluster_rows[0].keys() if cluster_rows else [])
                    if c != 'pdb_id']
    joined = []
    n_match = 0
    unmatched = []
    for r in summary_rows:
        key = r['pdb_id'].upper()
        match = by_key.get(key)
        merged = dict(r)
        if match:
            n_match += 1
            for c in cluster_cols:
                if c not in merged:
                    merged[c] = match[c]
        else:
            unmatched.append(r['pdb_id'])
            for c in cluster_cols:
                merged.setdefault(c, '')
        joined.append(merged)
    return joined, n_match, unmatched, cluster_cols


# ---------- Stats ----------

def ranks(values):
    """Average ranks (1-based) for a list of (index, value) pairs.

    Returns list aligned to original positions.
    """
    n = len(values)
    indexed = sorted(range(n), key=lambda i: values[i])
    out = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and values[indexed[j + 1]] == values[indexed[i]]:
            j += 1
        avg = (i + j) / 2 + 1  # 1-based average rank
        for k in range(i, j + 1):
            out[indexed[k]] = avg
        i = j + 1
    return out


def pearson(x, y):
    n = len(x)
    if n < 3:
        return None
    mx = sum(x) / n
    my = sum(y) / n
    num = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
    sx = math.sqrt(sum((xi - mx) ** 2 for xi in x))
    sy = math.sqrt(sum((yi - my) ** 2 for yi in y))
    if sx == 0 or sy == 0:
        return None
    return num / (sx * sy)


def spearman(x, y):
    """Returns (rho, p_two_sided) using normal approximation. None if undefined."""
    rho = pearson(ranks(x), ranks(y))
    if rho is None:
        return None, None
    n = len(x)
    if n < 4:
        return rho, None
    # t-statistic conversion → use normal approximation for p
    if abs(rho) >= 1.0:
        return rho, 0.0
    t = rho * math.sqrt((n - 2) / (1 - rho * rho))
    # Approximate p with two-sided normal; for n > ~30 the t and z agree well
    p = 2 * (1 - norm_cdf(abs(t)))
    return rho, p


def norm_cdf(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


def mann_whitney(group_a, group_b):
    """Two-sided Mann-Whitney U with normal approximation (tie correction)."""
    n1, n2 = len(group_a), len(group_b)
    if n1 == 0 or n2 == 0:
        return None, None, None
    combined = group_a + group_b
    r = ranks(combined)
    r1 = sum(r[:n1])
    u1 = r1 - n1 * (n1 + 1) / 2
    u2 = n1 * n2 - u1
    u = min(u1, u2)
    mu = n1 * n2 / 2
    # Tie correction
    from collections import Counter
    tied = Counter(combined)
    tie_term = sum(t * (t * t - 1) for t in tied.values())
    n = n1 + n2
    var = n1 * n2 / 12 * ((n + 1) - tie_term / (n * (n - 1))) if n > 1 else 0
    if var <= 0:
        return u, None, None
    sigma = math.sqrt(var)
    # continuity correction
    z = (u - mu + 0.5) / sigma if u < mu else (u - mu - 0.5) / sigma
    p = 2 * (1 - norm_cdf(abs(z)))
    return u, z, p


# ---------- Analysis ----------

CLUSTER_COLS = {'cluster_size', 'cluster_size_afdb50', 'cluster_size_full',
                'dedup_ratio', 'avg_pLDDT', 'avg_seq_len', 'nMem_file2'}


def collect_numeric_dim(joined, col):
    """Return aligned (f1_list, dim_list) for rows where both are present.

    For cluster-derived dimensions, drop rows with failed UniProt resolution
    (cluster_sizes.tsv 'NOT_FOUND' artifacts).
    """
    f1s, vals = [], []
    is_cluster_dim = col in CLUSTER_COLS
    for r in joined:
        f1 = as_float(r.get('f1'))
        v = as_float(r.get(col))
        if f1 is None or v is None:
            continue
        if is_cluster_dim and not cluster_lookup_ok(r):
            continue
        f1s.append(f1)
        vals.append(v)
    return f1s, vals


def categorical_groups(joined, col):
    """Group F1 by categorical value, dropping empty/NA cells."""
    from collections import defaultdict
    groups = defaultdict(list)
    for r in joined:
        f1 = as_float(r.get('f1'))
        cat = (r.get(col) or '').strip()
        if f1 is None or cat == '' or cat.upper() == 'NA':
            continue
        groups[cat].append(f1)
    return dict(groups)


def bucket_summary(joined, col, edges, labels):
    """For each bucket, return (label, mean F1, std F1, count)."""
    from collections import defaultdict
    buckets = defaultdict(list)
    is_cluster_dim = col in CLUSTER_COLS
    for r in joined:
        f1 = as_float(r.get('f1'))
        v = as_float(r.get(col))
        if f1 is None or v is None:
            continue
        if is_cluster_dim and not cluster_lookup_ok(r):
            continue
        for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
            if lo <= v <= hi if i == 0 else lo < v <= hi:
                buckets[labels[i]].append(f1)
                break
    rows = []
    for lab in labels:
        vals = buckets.get(lab, [])
        if not vals:
            rows.append((lab, None, None, 0))
            continue
        m = sum(vals) / len(vals)
        s = math.sqrt(sum((v - m) ** 2 for v in vals) / len(vals)) if len(vals) > 1 else 0
        rows.append((lab, m, s, len(vals)))
    return rows


def ascii_histogram(values, n_bins=20, width=55):
    if not values:
        return ['(empty)']
    lo, hi = min(values), max(values)
    if lo == hi:
        return [f'all values == {lo:.3f}']
    bw = (hi - lo) / n_bins
    counts = [0] * n_bins
    for v in values:
        idx = min(int((v - lo) / bw), n_bins - 1)
        counts[idx] += 1
    mc = max(counts)
    lines = []
    for i, c in enumerate(counts):
        edge = lo + i * bw
        bar = '█' * round(width * c / mc) if mc else ''
        lines.append(f'  {edge:>5.2f} |{bar} {c}')
    return lines


# ---------- Main ----------

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument('summary_tsv')
    ap.add_argument('--cluster-tsv', default='analysis/cluster_sizes.tsv')
    ap.add_argument('--out-dir', default='analysis')
    args = ap.parse_args()

    if not Path(args.summary_tsv).exists():
        sys.exit(f"ERROR: {args.summary_tsv} not found")
    if not Path(args.cluster_tsv).exists():
        sys.exit(f"ERROR: {args.cluster_tsv} not found")

    summary_rows = read_tsv(args.summary_tsv)
    cluster_rows = read_tsv(args.cluster_tsv)
    joined, n_match, unmatched, cluster_cols = join(summary_rows, cluster_rows)

    print("=" * 70)
    print(f"Inputs:  {args.summary_tsv} ({len(summary_rows)} rows)")
    print(f"         {args.cluster_tsv} ({len(cluster_rows)} rows)")
    print("=" * 70)
    print(f"Joined {n_match}/{len(joined)} PDBs; {len(unmatched)} unmatched")
    if unmatched:
        print(f"  Unmatched: {unmatched[:10]}{' …' if len(unmatched) > 10 else ''}")

    # Which expected dims are present?
    have_numeric = [c for c in NUMERIC_CANDIDATES if c in (summary_rows[0].keys() if summary_rows else []) or c in cluster_cols]
    have_categorical = [c for c in CATEGORICAL_CANDIDATES if c in cluster_cols]
    missing = [c for c in NUMERIC_CANDIDATES + CATEGORICAL_CANDIDATES
               if c not in set(have_numeric) | set(have_categorical)]
    if missing:
        print(f"\nMissing dimensions (re-run analyze_cluster_sizes.py with AFDB file 2 to enable): {missing}")

    # Cluster-lookup quality (cluster_sizes.tsv often has NOT_FOUND rows where the
    # bulk UniProt API failed; those need to be excluded from cluster-size stats).
    n_lookup_ok = sum(1 for r in joined if cluster_lookup_ok(r))
    n_lookup_failed = len(joined) - n_lookup_ok
    if n_lookup_failed:
        print(f"\nCluster lookup quality: {n_lookup_ok}/{len(joined)} PDBs have real cluster data; "
              f"{n_lookup_failed} have uniprot='NOT_FOUND' or cluster_size=0 and will be skipped "
              f"for cluster-derived correlations.")
        if n_lookup_ok < 100:
            print("  WARNING: too few clean cluster rows to draw reliable conclusions about cluster-size effects.")
        print("  → Re-run analyze_cluster_sizes.py — UniProt's bulk ID-mapping API likely rate-limited or timed out.")

    # Failure-mode classification: separate pipeline-failed-to-predict from predicted-wrong cohorts.
    print("\n" + "=" * 70)
    print("Failure-mode cohorts")
    print("=" * 70)
    no_pred = []         # n_predicted == 0 (alignment_mapping / pipeline failure)
    pred_zero_f1 = []    # n_predicted > 0 but F1 == 0 (predicted, all wrong)
    pred_partial = []    # 0 < F1 < 1
    pred_perfect = []    # F1 == 1
    for r in joined:
        f1 = as_float(r.get('f1'))
        np_ = as_float(r.get('n_predicted'))
        if f1 is None or np_ is None:
            continue
        if np_ == 0:
            no_pred.append(r)
        elif f1 == 0:
            pred_zero_f1.append(r)
        elif f1 >= 1.0:
            pred_perfect.append(r)
        else:
            pred_partial.append(r)
    total = len(no_pred) + len(pred_zero_f1) + len(pred_partial) + len(pred_perfect)
    print(f"  n_predicted=0 (no predictions at all):  {len(no_pred):>4}  ({100*len(no_pred)/total:.1f}%)")
    print(f"  predicted but F1=0 (all wrong):          {len(pred_zero_f1):>4}  ({100*len(pred_zero_f1)/total:.1f}%)")
    print(f"  partial success (0 < F1 < 1):            {len(pred_partial):>4}  ({100*len(pred_partial)/total:.1f}%)")
    print(f"  perfect (F1 = 1):                        {len(pred_perfect):>4}  ({100*len(pred_perfect)/total:.1f}%)")
    if no_pred:
        sample = [r['pdb_id'] for r in no_pred[:15]]
        print(f"\n  no-prediction PDBs (likely alignment_mapping or filter failures):")
        print(f"    {sample}{' …' if len(no_pred) > 15 else ''}")

    # Write joined wide table
    out_dir = Path(args.out_dir)
    out_table = out_dir / 'diagnostic_table.tsv'
    fieldnames = list(summary_rows[0].keys()) + [c for c in cluster_cols if c not in summary_rows[0]]
    write_tsv(out_table, joined, fieldnames)
    print(f"\nWrote {out_table}")

    # F1 distribution: ASCII histogram (visualize bimodality directly)
    f1_vals = [as_float(r['f1']) for r in joined]
    f1_vals = [v for v in f1_vals if v is not None]
    print("\n" + "=" * 70)
    print(f"F1 distribution  (n={len(f1_vals)}, mean={sum(f1_vals)/len(f1_vals):.3f}, "
          f"median={sorted(f1_vals)[len(f1_vals)//2]:.3f}, "
          f"zero F1: {sum(1 for v in f1_vals if v == 0)}, "
          f"perfect: {sum(1 for v in f1_vals if v == 1)})")
    print("=" * 70)
    for line in ascii_histogram(f1_vals):
        print(line)

    # Spearman ρ for each numeric dim
    print("\n" + "=" * 70)
    print("Spearman ρ:  f1 vs each numeric dimension (sorted by |ρ|)")
    print("=" * 70)
    spearman_rows = []
    for col in have_numeric:
        f1s, vals = collect_numeric_dim(joined, col)
        if len(f1s) < 5:
            continue
        rho, p = spearman(f1s, vals)
        if rho is None:
            continue
        spearman_rows.append((col, rho, p, len(f1s)))
    spearman_rows.sort(key=lambda r: -abs(r[1]))
    print(f"{'dimension':<25}{'ρ':>10}{'p (approx)':>14}{'n':>8}")
    for col, rho, p, n in spearman_rows:
        p_str = f'{p:.3g}' if p is not None else '—'
        print(f"{col:<25}{rho:>+10.3f}{p_str:>14}{n:>8}")

    # Mann-Whitney for categorical
    if have_categorical:
        print("\n" + "=" * 70)
        print("Mann-Whitney U:  f1 by categorical group")
        print("=" * 70)
        for col in have_categorical:
            groups = categorical_groups(joined, col)
            if len(groups) != 2:
                print(f"  {col}: {len(groups)} groups (need 2 for two-sample test). "
                      f"counts={[(k, len(v)) for k, v in groups.items()]}")
                continue
            keys = list(groups.keys())
            u, z, p = mann_whitney(groups[keys[0]], groups[keys[1]])
            mean_a = sum(groups[keys[0]]) / len(groups[keys[0]])
            mean_b = sum(groups[keys[1]]) / len(groups[keys[1]])
            p_str = f'{p:.3g}' if p is not None else '—'
            print(f"  {col}: {keys[0]} (n={len(groups[keys[0]])}, mean F1={mean_a:.3f}) vs "
                  f"{keys[1]} (n={len(groups[keys[1]])}, mean F1={mean_b:.3f})  "
                  f"U={u:.0f}  p={p_str}")

    # Bucket means
    print("\n" + "=" * 70)
    print("Bucket means")
    print("=" * 70)
    primary_cluster_col = next(
        (c for c in ('cluster_size_afdb50', 'cluster_size', 'cluster_size_full')
         if c in cluster_cols), None
    )
    if primary_cluster_col:
        edges = [-1, 20, 100, 500, float('inf')]
        labels = ['<20', '20-100', '100-500', '500+']
        bm = bucket_summary(joined, primary_cluster_col, edges, labels)
        print(f"\nF1 by {primary_cluster_col}:")
        print(f"  {'bucket':<12}{'mean':>10}{'std':>10}{'n':>8}")
        for lab, mean, std, n in bm:
            mean_s = f'{mean:.3f}' if mean is not None else '—'
            std_s = f'{std:.3f}' if std is not None else '—'
            print(f"  {lab:<12}{mean_s:>10}{std_s:>10}{n:>8}")

    if 'n_true' in (summary_rows[0].keys() if summary_rows else []):
        edges = [0, 1, 3, 5, float('inf')]
        labels = ['1', '2-3', '4-5', '6+']
        bm = bucket_summary(joined, 'n_true', edges, labels)
        print("\nF1 by n_true (M-CSA ground-truth count):")
        print(f"  {'bucket':<12}{'mean':>10}{'std':>10}{'n':>8}")
        for lab, mean, std, n in bm:
            mean_s = f'{mean:.3f}' if mean is not None else '—'
            std_s = f'{std:.3f}' if std is not None else '—'
            print(f"  {lab:<12}{mean_s:>10}{std_s:>10}{n:>8}")

    # Long-format TSV for easy plotting in any tool
    long_path = out_dir / 'diagnostic_long.tsv'
    Path(long_path).parent.mkdir(parents=True, exist_ok=True)
    with open(long_path, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['pdb_id', 'dimension', 'value', 'f1'])
        for col in have_numeric:
            for r in joined:
                f1 = as_float(r.get('f1'))
                v = as_float(r.get(col))
                if f1 is None or v is None:
                    continue
                w.writerow([r['pdb_id'], col, v, f1])
    print(f"\nWrote {long_path}  (long-format for plotting)")

    # Recommendation
    print("\n" + "=" * 70)
    if spearman_rows:
        top_col, top_rho, top_p, top_n = spearman_rows[0]
        moderate = [r for r in spearman_rows if 0.1 <= abs(r[1]) <= 0.3]
        if abs(top_rho) > 0.3:
            print(f"RECOMMENDATION (Branch 1 — tunable): {top_col} shows |ρ|={abs(top_rho):.2f}.")
            print("  Failure mode partially explained by this dimension; consider filtering")
            print("  the benchmark to a usable regime or conditioning scoring on it.")
        elif len(moderate) >= 2:
            cols = ', '.join(c for c, *_ in moderate)
            print(f"RECOMMENDATION (Branch 2 — multi-factorial): {len(moderate)} dimensions "
                  f"with 0.1 ≤ |ρ| ≤ 0.3 ({cols}).")
            print("  Failure is multi-factorial — next step is multivariate regression or")
            print("  per-residue rank diagnostics (extend analyze_top_n_comparison.py).")
        else:
            print(f"RECOMMENDATION (Branch 3 — structural): no cluster-property correlation "
                  f"above |ρ|={abs(top_rho):.2f}.")
            print("  Bimodality is not explained by cluster size/quality. Extend")
            print("  analyze_top_n_comparison.py to study residue-level ranking failures —")
            print("  likely an alignment-mapping or scoring-formula issue, not cluster quality.")
    else:
        print("No correlations could be computed — check the join and that numeric columns parse.")
    print("=" * 70)


if __name__ == '__main__':
    main()
