"""
Parse foldseek easy-search TSV output into a UniProt accession list.

Expects the TSV to have been produced with:
    --format-output query,target,fident,alnlen,evalue,prob,lddt,alntmscore

Filters by minimum probability, takes top-N by probability (descending),
extracts the UniProt accession from each target ID
(AFDB targets look like 'AF-P12345-F1-model_v4'), and writes one
accession per line to stdout or to -o.

Exits non-zero if no hits pass the filter — preferring a hard fail over
silently producing an empty MSA downstream.

Usage:
    python3 src/parse_foldseek_hits.py <foldseek_hits.tsv>
        [--prob-min 0.8] [--top-n 200] [-o homologs.txt]
"""

import argparse
import csv
import re
import sys
from pathlib import Path


AF_TARGET_RE = re.compile(r"AF-([A-Z0-9]+)-F\d+-model_v\d+")

COLUMNS = ["query", "target", "fident", "alnlen", "evalue",
           "prob", "lddt", "alntmscore"]


def parse_hits(tsv_path, prob_min, top_n):
    """Return a list of (uniprot, prob) tuples, top-N by prob desc."""
    kept = []
    n_total = 0
    n_low_prob = 0
    n_unparsable = 0

    with open(tsv_path) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            n_total += 1
            if len(row) < len(COLUMNS):
                n_unparsable += 1
                continue
            target = row[1]
            try:
                prob = float(row[5])
            except ValueError:
                n_unparsable += 1
                continue

            if prob < prob_min:
                n_low_prob += 1
                continue

            m = AF_TARGET_RE.search(target)
            if not m:
                n_unparsable += 1
                continue

            kept.append((m.group(1), prob))

    kept.sort(key=lambda x: x[1], reverse=True)
    truncated = kept[:top_n]

    print(f"  Foldseek hits: {n_total} total, "
          f"{n_low_prob} below prob {prob_min}, "
          f"{n_unparsable} unparsable, "
          f"{len(kept)} kept (capped at {top_n} → {len(truncated)})",
          file=sys.stderr)

    return truncated


def main():
    p = argparse.ArgumentParser(description="Parse foldseek hits → UniProt list")
    p.add_argument("tsv", type=Path, help="foldseek easy-search TSV output")
    p.add_argument("--prob-min", type=float, default=0.8,
                   help="Minimum prob to keep a hit (default: 0.8)")
    p.add_argument("--top-n", type=int, default=200,
                   help="Maximum hits to keep, by prob descending (default: 200)")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="Output file (default: stdout)")
    args = p.parse_args()

    if not args.tsv.exists():
        print(f"ERROR: foldseek TSV not found at {args.tsv}", file=sys.stderr)
        sys.exit(1)

    hits = parse_hits(args.tsv, args.prob_min, args.top_n)

    if not hits:
        print(f"ERROR: no foldseek hits passed prob >= {args.prob_min}. "
              f"Cannot build MSA from an empty homolog list.", file=sys.stderr)
        sys.exit(1)

    out_lines = [u for u, _ in hits]
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text("\n".join(out_lines) + "\n")
        print(f"  Wrote {len(out_lines)} UniProt accessions → {args.output}",
              file=sys.stderr)
    else:
        sys.stdout.write("\n".join(out_lines) + "\n")


if __name__ == "__main__":
    main()
