#!/usr/bin/env python3
"""
Extract a PDB-id list from a parsed M-CSA TSV.

Reads is_reference=True entries and writes a TSV:
  mcsa_id  pdb_id  n_catalytic_residues

Usage:
    python3 extract_pdb_list_from_mcsa.py <catalytic_residues_homologues_parsed.tsv> \
        [-o mcsa_representatives_parsed_monomers.tsv] [--max N] [--shuffle]
"""

import csv
import argparse
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(description='Extract a PDB-id list from a parsed M-CSA TSV')
    parser.add_argument('mcsa_tsv', help='Parsed M-CSA TSV file')
    parser.add_argument('-o', '--output', default='mcsa_representatives_parsed_monomers.tsv',
                        help='Output TSV file')
    parser.add_argument('--max', type=int, default=None, help='Max entries (for subset runs)')
    parser.add_argument('--min-residues', type=int, default=1, help='Min catalytic residues to include')
    parser.add_argument('--shuffle', action='store_true', help='Randomize order (for subset sampling)')
    
    args = parser.parse_args()
    
    # Collect all reference entries
    # Group by mcsa_id to handle multi-chain (take first chain)
    entries = defaultdict(lambda: {'pdb_id': None, 'residues': set()})
    
    with open(args.mcsa_tsv, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        # Detect column name
        fieldnames = reader.fieldnames
        res_col = 'chain_residues' if 'chain_residues' in fieldnames else 'residues'
        
        for row in reader:
            if row['is_reference'].strip().upper() != 'TRUE':
                continue
            
            mcsa_id = row['mcsa_id'].strip()
            pdb_id = row['pdb_id'].strip().upper()
            residues_str = row[res_col].strip()
            
            entries[mcsa_id]['pdb_id'] = pdb_id
            
            # Count non-gap residues
            for token in residues_str.split(','):
                token = token.strip()
                if token and token != '-':
                    entries[mcsa_id]['residues'].add(token)
    
    # Build output list
    results = []
    for mcsa_id, info in sorted(entries.items(), key=lambda x: int(x[0]) if x[0].isdigit() else 0):
        n_res = len(info['residues'])
        if n_res < args.min_residues:
            continue
        results.append((mcsa_id, info['pdb_id'], n_res))
    
    if args.shuffle:
        import random
        random.shuffle(results)
    
    if args.max:
        results = results[:args.max]
    
    # Write output
    with open(args.output, 'w') as f:
        f.write("mcsa_id\tpdb_id\tn_catalytic_residues\n")
        for mcsa_id, pdb_id, n_res in results:
            f.write(f"{mcsa_id}\t{pdb_id}\t{n_res}\n")
    
    print(f"Generated benchmark list: {len(results)} entries")
    print(f"Output: {args.output}")
    
    # Stats
    n_res_list = [r[2] for r in results]
    if n_res_list:
        print(f"Catalytic residues per entry: min={min(n_res_list)}, max={max(n_res_list)}, "
              f"mean={sum(n_res_list)/len(n_res_list):.1f}")


if __name__ == '__main__':
    main()
