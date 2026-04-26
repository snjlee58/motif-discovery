#!/usr/bin/env python3
"""
Analyze AFDB cluster sizes for the PDBs in a benchmark TSV.

For each PDB:
  1. Resolve PDB → UniProt via UniProt ID Mapping API (one bulk call).
  2. Look up the UniProt's AFDB cluster representative (file 5).
  3. Count cluster members, broken down by cluFlag:
       - cluFlag=1: AFDB50 reps (Foldseek-clustered) — the 50%-identity dedup
       - cluFlag=2: non-rep AFDB50 members (near-duplicates of their AFDB50 rep)
  4. Optionally enrich with file 2 metadata (avg pLDDT, LCA taxon, etc.)

Outputs:
  - <output>.tsv with per-PDB cluster sizes + dedup ratio + metadata
  - Summary stats + ASCII histograms (full vs AFDB50-rep-only)

Usage (from repo root):
    python3 analysis/analyze_cluster_sizes.py \
        benchmark/mcsa_representatives_parsed_monomers.tsv \
        --output analysis/cluster_sizes.tsv

Cluster files default to $FAST/afdb_clusters/v6/. If file 2 (cluster overview) is
missing, metadata columns are skipped — only the cluster-size analysis runs.
File 2 download (v6 — match AlphaFold v6 structures used by pipeline.sh):
    cd $FAST/afdb_clusters/v6
    wget https://afdb-cluster.steineggerlab.workers.dev/v6/2-repId_isDark_nMem_repLen_avgLen_repPlddt_avgPlddt_LCAtaxId.tsv.gz
    gunzip 2-repId_isDark_nMem_repLen_avgLen_repPlddt_avgPlddt_LCAtaxId.tsv.gz
"""

import argparse
import csv
import json
import os
import re
import sys
import time
from collections import Counter
from pathlib import Path
from urllib import parse, request


UNIPROT_RUN = "https://rest.uniprot.org/idmapping/run"
UNIPROT_STATUS = "https://rest.uniprot.org/idmapping/status/"
UNIPROT_RESULTS = "https://rest.uniprot.org/idmapping/results/"
PAGE_SIZE = 500  # UniProt's max page size for paginated results


def submit_pdb_to_uniprot(pdb_ids):
    """Submit a batch PDB → UniProtKB mapping job. Returns job ID."""
    data = parse.urlencode({
        'from': 'PDB',
        'to': 'UniProtKB',
        'ids': ','.join(pdb_ids),
    }).encode()
    req = request.Request(UNIPROT_RUN, data=data)
    with request.urlopen(req, timeout=30) as resp:
        return json.loads(resp.read())['jobId']


def parse_results(results):
    """Parse UniProt API results into {pdb_lower: uniprot_acc}."""
    mapping = {}
    for entry in results:
        pdb = entry['from'].lower()
        to = entry['to']
        if isinstance(to, dict):
            uniprot = to.get('primaryAccession', '')
        else:
            uniprot = to
        if uniprot and pdb not in mapping:  # take first if multiple
            mapping[pdb] = uniprot
    return mapping


def fetch_results(job_id, max_wait=180):
    """Wait for job, then fetch ALL results via paginated endpoint.

    /idmapping/status/ silently caps responses at 25 entries (default page size),
    so we wait for completion and then walk the paginated /results/ endpoint with
    size=500, following the Link header for additional pages.
    """
    start = time.time()
    while time.time() - start < max_wait:
        try:
            with request.urlopen(f"{UNIPROT_STATUS}{job_id}", timeout=30) as resp:
                data = json.loads(resp.read())
            status = data.get('jobStatus')
            if status in (None, 'FINISHED') or 'results' in data or 'redirectURL' in data:
                break
            if status in ('RUNNING', 'NEW'):
                time.sleep(2)
                continue
        except Exception as e:
            print(f"  status check failed: {e}", file=sys.stderr)
            time.sleep(2)
    else:
        raise TimeoutError(f"UniProt job {job_id} did not complete in {max_wait}s")

    all_results = []
    url = f"{UNIPROT_RESULTS}{job_id}?format=json&size={PAGE_SIZE}"
    page = 0
    while url:
        page += 1
        with request.urlopen(url, timeout=60) as resp:
            data = json.loads(resp.read())
            link_header = resp.headers.get('Link', '')
        all_results.extend(data.get('results', []))
        m = re.search(r'<([^>]+)>;\s*rel="next"', link_header)
        url = m.group(1) if m else None

    print(f"  fetched {len(all_results)} mappings ({page} page{'s' if page != 1 else ''})",
          file=sys.stderr)
    return parse_results(all_results)


def scan_cluster_file(cluster_file, target_uniprots):
    """Single pass through file 5.

    Counts members per rep, split by cluFlag:
      - cluFlag=1: AFDB50 representatives (Foldseek-clustered) — the 50%-identity dedup
      - cluFlag=2: non-rep AFDB50 members (near-duplicates of an AFDB50 rep)

    Returns:
      sizes_full: Counter {rep_id: total member count}
      sizes_afdb50: Counter {rep_id: cluFlag=1 count}
      uniprot_to_rep: {uniprot: rep_id} for targets only
    """
    sizes_full = Counter()
    sizes_afdb50 = Counter()
    uniprot_to_rep = {}

    print(f"Scanning {cluster_file}...", file=sys.stderr)
    n_rows = 0
    with open(cluster_file) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            rep, mem, clu_flag = parts[0], parts[1], parts[2]
            sizes_full[rep] += 1
            if clu_flag == '1':
                sizes_afdb50[rep] += 1
            if mem in target_uniprots and mem not in uniprot_to_rep:
                uniprot_to_rep[mem] = rep
            n_rows += 1
            if n_rows % 20_000_000 == 0:
                print(f"  {n_rows:>12,} rows · {len(uniprot_to_rep)}/{len(target_uniprots)} targets found",
                      file=sys.stderr)

    print(f"Done. {n_rows:,} rows · {len(sizes_full):,} unique clusters.",
          file=sys.stderr)
    return sizes_full, sizes_afdb50, uniprot_to_rep


def load_cluster_overview(file_path, target_reps):
    """Load file 2 cluster overview for target reps.

    Format: repId, isDark, nMem, repLen, avgLen, repPlddt, avgPlddt, LCAtaxID

    Returns: {rep_id: {nMem, isDark, avgLen, avgPlddt, LCAtaxID}} or {} if file missing.
    """
    if not Path(file_path).exists():
        print(f"File 2 not found at {file_path} — skipping metadata enrichment",
              file=sys.stderr)
        return {}

    print(f"Loading cluster metadata from {file_path}...", file=sys.stderr)
    overview = {}
    with open(file_path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 8:
                continue
            rep = parts[0]
            if rep in target_reps:
                try:
                    overview[rep] = {
                        'isDark': parts[1],
                        'nMem': int(parts[2]),
                        'avgLen': float(parts[4]),
                        'avgPlddt': float(parts[6]),
                        'LCAtaxID': parts[7],
                    }
                except (ValueError, IndexError):
                    continue
    print(f"  metadata for {len(overview)}/{len(target_reps)} target clusters",
          file=sys.stderr)
    return overview


def histogram(values):
    """Print an ASCII histogram of cluster sizes (AFDB-paper-style bins)."""
    bins = [(2, 4), (5, 10), (11, 20), (21, 40), (41, 80),
            (81, 160), (161, 320), (321, 700), (701, 1800),
            (1801, 10000), (10001, 100000)]

    counts = [(lo, hi, sum(1 for v in values if lo <= v <= hi)) for lo, hi in bins]
    max_count = max((c for _, _, c in counts), default=1) or 1
    bar_max = 40

    print()
    print(f"{'Cluster size':<20}{'Count':>8}  Distribution")
    print("=" * 70)
    for lo, hi, n in counts:
        bar = '█' * round(bar_max * n / max_count)
        print(f"  {lo:>5}-{hi:<10}{n:>8}  {bar}")


def stats(label, values):
    """Print summary stats for a list of cluster sizes."""
    if not values:
        return
    s = sorted(values)
    n = len(s)
    print(f"{label}:")
    print(f"  Min:    {min(s):>10,}")
    print(f"  Median: {s[n // 2]:>10,}")
    print(f"  Mean:   {sum(s) / n:>10,.1f}")
    print(f"  Max:    {max(s):>10,}")
    print(f"  P90:    {s[int(n * 0.9)]:>10,}")
    print(f"  P95:    {s[int(n * 0.95)]:>10,}")
    print(f"  P99:    {s[int(n * 0.99)]:>10,}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('tsv', help='Benchmark TSV (mcsa_id, pdb_id, n_catalytic_residues)')
    parser.add_argument('--cluster-file',
                        default=os.path.join(
                            os.environ.get('FAST', '/fast/sunny'),
                            'afdb_clusters/v6/5-allmembers-repId-entryId-cluFlag-taxId.tsv'),
                        help='AFDB file 5: all members')
    parser.add_argument('--cluster-overview',
                        default=os.path.join(
                            os.environ.get('FAST', '/fast/sunny'),
                            'afdb_clusters/v6/2-repId_isDark_nMem_repLen_avgLen_repPlddt_avgPlddt_LCAtaxId.tsv'),
                        help='AFDB file 2: per-cluster overview (optional)')
    parser.add_argument('--output', default='analysis/cluster_sizes.tsv',
                        help='Output TSV (default: analysis/cluster_sizes.tsv)')
    parser.add_argument('--top', type=int, default=20,
                        help='Print top N largest clusters in summary')
    args = parser.parse_args()

    if not Path(args.cluster_file).exists():
        sys.exit(f"ERROR: cluster file not found: {args.cluster_file}")

    # Read benchmark TSV
    proteins = []
    with open(args.tsv) as f:
        next(f)
        for line in f:
            mcsa_id, pdb, n_res = line.strip().split('\t')
            proteins.append({'mcsa_id': mcsa_id, 'pdb': pdb.upper(), 'n_res': int(n_res)})
    print(f"Loaded {len(proteins)} PDBs from {args.tsv}", file=sys.stderr)

    # Bulk PDB → UniProt
    print("Resolving PDB → UniProt (bulk)...", file=sys.stderr)
    job_id = submit_pdb_to_uniprot([p['pdb'].lower() for p in proteins])
    print(f"  job: {job_id}", file=sys.stderr)
    pdb_to_uniprot = fetch_results(job_id)
    print(f"  resolved {len(pdb_to_uniprot)}/{len(proteins)} PDBs", file=sys.stderr)

    target_uniprots = set(pdb_to_uniprot.values())

    # Single-pass file 5 scan with cluFlag breakdown
    sizes_full, sizes_afdb50, uniprot_to_rep = scan_cluster_file(
        args.cluster_file, target_uniprots)

    # File 2 metadata (optional — gracefully skipped if absent)
    target_reps = set(uniprot_to_rep.values())
    overview = load_cluster_overview(args.cluster_overview, target_reps)

    # Build output rows
    rows = []
    full_vals = []
    afdb50_vals = []
    not_in_afdb = []
    no_uniprot = []
    for p in proteins:
        pdb = p['pdb']
        uniprot = pdb_to_uniprot.get(pdb.lower(), '')
        if not uniprot:
            no_uniprot.append(pdb)
            rep = ''
            full_size = afdb50_size = 0
        else:
            rep = uniprot_to_rep.get(uniprot, '')
            full_size = sizes_full.get(rep, 0) if rep else 0
            afdb50_size = sizes_afdb50.get(rep, 0) if rep else 0
            if full_size == 0:
                not_in_afdb.append(pdb)

        meta = overview.get(rep, {})
        ratio = (full_size / afdb50_size) if afdb50_size else 0.0
        rows.append({
            'mcsa_id': p['mcsa_id'],
            'pdb_id': pdb,
            'uniprot': uniprot or 'NOT_FOUND',
            'rep_id': rep or 'NOT_FOUND',
            'cluster_size_full': full_size,
            'cluster_size_afdb50': afdb50_size,
            'dedup_ratio': f"{ratio:.2f}" if ratio else '',
            'nMem_file2': meta.get('nMem', ''),
            'avg_seq_len': f"{meta.get('avgLen'):.1f}" if 'avgLen' in meta else '',
            'avg_pLDDT': f"{meta.get('avgPlddt'):.1f}" if 'avgPlddt' in meta else '',
            'isDark': meta.get('isDark', ''),
            'LCA_taxID': meta.get('LCAtaxID', ''),
        })
        if full_size > 0:
            full_vals.append(full_size)
            if afdb50_size > 0:
                afdb50_vals.append(afdb50_size)

    rows.sort(key=lambda r: -r['cluster_size_full'])

    # Write output
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ['mcsa_id', 'pdb_id', 'uniprot', 'rep_id',
                  'cluster_size_full', 'cluster_size_afdb50', 'dedup_ratio',
                  'nMem_file2', 'avg_seq_len', 'avg_pLDDT', 'isDark', 'LCA_taxID']
    with open(args.output, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {args.output}", file=sys.stderr)

    # Summary
    if full_vals:
        n = len(full_vals)
        print()
        print("=" * 70)
        print(f"Cluster size summary  ({n}/{len(proteins)} PDBs resolved to AFDB cluster)")
        print("=" * 70)
        stats("Full cluster size (file 5, all members)", full_vals)
        print()
        stats("AFDB50-rep-only size (cluFlag=1, the 'free dedup')", afdb50_vals)

        # Dedup ratios
        ratios = [a / b for a, b in zip(full_vals, afdb50_vals) if b > 0]
        if ratios:
            rs = sorted(ratios)
            n = len(rs)
            print()
            print(f"Dedup ratio (full / afdb50) — how much you save by filtering to cluFlag=1:")
            print(f"  Median: {rs[n // 2]:>6.2f}x")
            print(f"  Mean:   {sum(rs) / n:>6.2f}x")
            print(f"  Max:    {max(rs):>6.2f}x")

        # Histograms — both for direct comparison
        print()
        print("Distribution: FULL cluster size")
        histogram(full_vals)
        print()
        print("Distribution: AFDB50-rep-only size")
        histogram(afdb50_vals)

        print()
        print(f"Top {args.top} largest (sorted by full size):")
        print(f"  {'PDB':<8}{'UniProt':<12}{'Full':>10}{'AFDB50':>10}{'Ratio':>8}{'avgPlddt':>10}")
        for r in rows[:args.top]:
            ratio_str = r['dedup_ratio'] or '-'
            plddt = r['avg_pLDDT'] or '-'
            print(f"  {r['pdb_id']:<8}{r['uniprot']:<12}"
                  f"{r['cluster_size_full']:>10,}{r['cluster_size_afdb50']:>10,}"
                  f"{ratio_str:>8}{plddt:>10}")

    if no_uniprot:
        print(f"\n{len(no_uniprot)} PDBs failed UniProt resolution: "
              f"{no_uniprot[:10]}{'...' if len(no_uniprot) > 10 else ''}",
              file=sys.stderr)
    if not_in_afdb:
        print(f"{len(not_in_afdb)} PDBs resolved but not in AFDB clusters: "
              f"{not_in_afdb[:10]}{'...' if len(not_in_afdb) > 10 else ''}",
              file=sys.stderr)


if __name__ == '__main__':
    main()
