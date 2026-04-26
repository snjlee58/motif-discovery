#!/usr/bin/env python3
"""
Analyze AFDB cluster sizes for the PDBs in a benchmark TSV.

For each PDB:
  1. Resolve PDB → UniProt via UniProt ID Mapping API (one bulk call).
  2. Look up the UniProt's AFDB cluster representative.
  3. Count cluster members.

Outputs:
  - <output>.tsv with per-PDB cluster sizes (sorted largest first).
  - Summary stats + ASCII histogram printed to stdout.

Usage (from repo root):
    python3 analysis/analyze_cluster_sizes.py \
        benchmark/mcsa_representatives_parsed_monomers.tsv \
        --output analysis/cluster_sizes.tsv

The cluster file path defaults to $FAST/afdb_clusters/5-allmembers-...tsv
(same convention as pipeline.sh).
"""

import argparse
import csv
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path
from urllib import parse, request


UNIPROT_RUN = "https://rest.uniprot.org/idmapping/run"
UNIPROT_STATUS = "https://rest.uniprot.org/idmapping/status/"
UNIPROT_RESULTS = "https://rest.uniprot.org/idmapping/results/"


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
    """Poll for job completion, fetch results. Returns {pdb_lower: uniprot}."""
    start = time.time()
    while time.time() - start < max_wait:
        try:
            with request.urlopen(f"{UNIPROT_STATUS}{job_id}", timeout=30) as resp:
                data = json.loads(resp.read())
            if data.get('results'):
                return parse_results(data['results'])
            if data.get('jobStatus') == 'RUNNING':
                time.sleep(2)
                continue
        except Exception as e:
            print(f"  status check failed: {e}", file=sys.stderr)

        # Fallback: stream endpoint
        try:
            url = f"{UNIPROT_RESULTS}stream/{job_id}?format=json"
            with request.urlopen(url, timeout=60) as resp:
                data = json.loads(resp.read())
            if data.get('results'):
                return parse_results(data['results'])
        except Exception as e:
            print(f"  stream fetch failed: {e}", file=sys.stderr)

        time.sleep(2)

    raise TimeoutError(f"UniProt job {job_id} did not complete in {max_wait}s")


def scan_cluster_file(cluster_file, target_uniprots):
    """Single pass through the cluster file.

    Returns:
      cluster_sizes: Counter {rep_id: member_count} for ALL clusters
      uniprot_to_rep: {uniprot: rep_id} for target uniprots only
    """
    cluster_sizes = Counter()
    uniprot_to_rep = {}

    print(f"Scanning {cluster_file}...", file=sys.stderr)
    n_rows = 0
    with open(cluster_file) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            rep, mem = parts[0], parts[1]
            cluster_sizes[rep] += 1
            if mem in target_uniprots and mem not in uniprot_to_rep:
                uniprot_to_rep[mem] = rep
            n_rows += 1
            if n_rows % 20_000_000 == 0:
                print(f"  {n_rows:>12,} rows · {len(uniprot_to_rep)}/{len(target_uniprots)} targets found",
                      file=sys.stderr)

    print(f"Done. {n_rows:,} rows · {len(cluster_sizes):,} unique clusters.",
          file=sys.stderr)
    return cluster_sizes, uniprot_to_rep


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('tsv', help='Benchmark TSV (mcsa_id, pdb_id, n_catalytic_residues)')
    parser.add_argument('--cluster-file',
                        default=os.path.join(
                            os.environ.get('FAST', '/fast/sunny'),
                            'afdb_clusters/5-allmembers-repId-entryId-cluFlag-taxId.tsv'),
                        help='AFDB 5-allmembers cluster file')
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

    # Single-pass cluster file scan
    cluster_sizes, uniprot_to_rep = scan_cluster_file(args.cluster_file, target_uniprots)

    # Build output rows
    rows = []
    sizes = []
    not_in_afdb = []
    no_uniprot = []
    for p in proteins:
        pdb = p['pdb']
        uniprot = pdb_to_uniprot.get(pdb.lower(), '')
        if not uniprot:
            no_uniprot.append(pdb)
            rep, size = '', 0
        else:
            rep = uniprot_to_rep.get(uniprot, '')
            size = cluster_sizes.get(rep, 0) if rep else 0
            if size == 0:
                not_in_afdb.append(pdb)
        rows.append({
            'mcsa_id': p['mcsa_id'],
            'pdb_id': pdb,
            'uniprot': uniprot or 'NOT_FOUND',
            'rep_id': rep or 'NOT_FOUND',
            'cluster_size': size,
        })
        if size > 0:
            sizes.append(size)

    rows.sort(key=lambda r: -r['cluster_size'])

    # Write output
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['mcsa_id', 'pdb_id', 'uniprot', 'rep_id', 'cluster_size'],
                           delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {args.output}", file=sys.stderr)

    # Summary
    if sizes:
        s = sorted(sizes)
        n = len(s)
        print()
        print("=" * 70)
        print(f"Cluster size summary  ({n}/{len(proteins)} PDBs resolved to AFDB cluster)")
        print("=" * 70)
        print(f"  Min:    {min(s):>10,}")
        print(f"  Median: {s[n // 2]:>10,}")
        print(f"  Mean:   {sum(s) / n:>10,.1f}")
        print(f"  Max:    {max(s):>10,}")
        print(f"  P90:    {s[int(n * 0.9)]:>10,}")
        print(f"  P95:    {s[int(n * 0.95)]:>10,}")
        print(f"  P99:    {s[int(n * 0.99)]:>10,}")
        histogram(sizes)

        print()
        print(f"Top {args.top} largest clusters:")
        print(f"  {'PDB':<8}{'UniProt':<12}{'Cluster size':>14}")
        for r in rows[:args.top]:
            print(f"  {r['pdb_id']:<8}{r['uniprot']:<12}{r['cluster_size']:>14,}")

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
