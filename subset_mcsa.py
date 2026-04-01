#!/usr/bin/env python3
"""
Analyze M-CSA reference entries and create subsets by complexity.

Reads the parsed TSV and categorizes entries by:
  - Number of chains (single vs multi-chain in ground truth)
  - Number of catalytic residues
  - Whether all residues are on one chain

Usage:
    python3 subset_mcsa.py <catalytic_residues_separated_aligned.tsv> [-o output_dir]
"""

import csv
import argparse
from collections import defaultdict
from pathlib import Path


def analyze_entry(rows):
    """Analyze a set of rows for one mcsa_id + pdb_id combination."""
    pdb_id = rows[0]['pdb_id'].strip().upper()
    mcsa_id = rows[0]['mcsa_id'].strip()
    
    res_col = 'chain_residues' if 'chain_residues' in rows[0] else 'residues'
    
    # Collect all residues and chains
    all_residues = set()
    chains_seen = set()
    residues_per_chain = defaultdict(set)
    
    for row in rows:
        residues_str = row[res_col].strip()
        for token in residues_str.split(','):
            token = token.strip()
            if token == '-' or not token:
                continue
            chain = token[0]
            try:
                resnum = int(token[1:])
            except ValueError:
                continue
            chains_seen.add(chain)
            all_residues.add(token)
            residues_per_chain[chain].add(resnum)
    
    n_chains_in_gt = len(chains_seen)
    n_rows = len(rows)  # number of chain rows for this pdb
    n_residues = len(all_residues)
    
    # Get unique residue numbers (ignoring chain letter)
    unique_resnums = set()
    for token in all_residues:
        try:
            unique_resnums.add(int(token[1:]))
        except ValueError:
            pass
    
    return {
        'mcsa_id': mcsa_id,
        'pdb_id': pdb_id,
        'n_chains_in_gt': n_chains_in_gt,
        'n_rows': n_rows,
        'n_residues_total': n_residues,
        'n_unique_resnums': len(unique_resnums),
        'chains': sorted(chains_seen),
        'is_single_chain': n_chains_in_gt == 1 and n_rows == 1,
    }


def main():
    parser = argparse.ArgumentParser(description='Analyze and subset M-CSA entries')
    parser.add_argument('mcsa_tsv', help='Parsed M-CSA TSV file')
    parser.add_argument('-o', '--output-dir', default='.', help='Output directory for subsets')
    
    args = parser.parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    res_col = None
    
    # Read all reference entries, grouped by (mcsa_id, pdb_id)
    entries_by_key = defaultdict(list)
    
    with open(args.mcsa_tsv, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        # Detect column name
        res_col = 'chain_residues' if 'chain_residues' in reader.fieldnames else 'residues'
        
        for row in reader:
            if row['is_reference'].strip().upper() != 'TRUE':
                continue
            key = (row['mcsa_id'].strip(), row['pdb_id'].strip().upper())
            entries_by_key[key].append(row)
    
    # Analyze each entry
    analyses = []
    for key, rows in entries_by_key.items():
        info = analyze_entry(rows)
        analyses.append(info)
    
    # Print overall stats
    total = len(analyses)
    single_chain = [a for a in analyses if a['is_single_chain']]
    multi_chain = [a for a in analyses if not a['is_single_chain']]
    
    print(f"{'='*60}")
    print(f"M-CSA Reference Entry Analysis")
    print(f"{'='*60}")
    print(f"Total reference entries: {total}")
    print(f"Single-chain (1 chain, 1 row): {len(single_chain)} ({len(single_chain)/total*100:.1f}%)")
    print(f"Multi-chain: {len(multi_chain)} ({len(multi_chain)/total*100:.1f}%)")
    
    # Residue count distribution
    res_counts = [a['n_unique_resnums'] for a in analyses]
    print(f"\nCatalytic residues per entry:")
    print(f"  Min: {min(res_counts)}, Max: {max(res_counts)}, "
          f"Mean: {sum(res_counts)/len(res_counts):.1f}, "
          f"Median: {sorted(res_counts)[len(res_counts)//2]}")
    
    # Buckets by residue count
    print(f"\nDistribution by catalytic residue count:")
    for lo, hi, label in [(1,3,'1-3 (small)'), (4,6,'4-6 (medium)'), 
                           (7,10,'7-10 (large)'), (11,999,'11+ (very large)')]:
        count = sum(1 for r in res_counts if lo <= r <= hi)
        print(f"  {label}: {count} ({count/total*100:.1f}%)")
    
    # Create subsets
    subsets = {
        'single_chain': single_chain,
        'single_chain_small': [a for a in single_chain if a['n_unique_resnums'] <= 6],
        'single_chain_medium': [a for a in single_chain if 4 <= a['n_unique_resnums'] <= 8],
        'single_chain_large': [a for a in single_chain if a['n_unique_resnums'] >= 7],
        'multi_chain': multi_chain,
        'all': analyses,
    }
    
    print(f"\n{'='*60}")
    print(f"Subsets created:")
    print(f"{'='*60}")
    
    for name, entries in subsets.items():
        if not entries:
            continue
        
        outfile = output_dir / f"benchmark_{name}.tsv"
        with open(outfile, 'w') as f:
            f.write("mcsa_id\tpdb_id\tn_catalytic_residues\n")
            for a in sorted(entries, key=lambda x: int(x['mcsa_id']) if x['mcsa_id'].isdigit() else 0):
                f.write(f"{a['mcsa_id']}\t{a['pdb_id']}\t{a['n_unique_resnums']}\n")
        
        res_counts_sub = [a['n_unique_resnums'] for a in entries]
        print(f"  {name}: {len(entries)} entries "
              f"(avg {sum(res_counts_sub)/len(res_counts_sub):.1f} cat residues) "
              f"→ {outfile}")
    
    # Print some examples from single_chain_small for quick testing
    print(f"\n{'='*60}")
    print(f"Sample entries for quick testing (single_chain_small):")
    print(f"{'='*60}")
    print(f"{'MCSA':<8} {'PDB':<8} {'#Res':<6} {'Chains'}")
    print(f"{'-'*35}")
    for a in sorted(subsets['single_chain_small'], 
                     key=lambda x: x['n_unique_resnums'])[:20]:
        print(f"{a['mcsa_id']:<8} {a['pdb_id']:<8} {a['n_unique_resnums']:<6} {','.join(a['chains'])}")


if __name__ == '__main__':
    main()