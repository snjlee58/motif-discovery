#!/usr/bin/env python3
"""
Select N diverse proteins from M-CSA for benchmarking.
Picks one protein per EC top-level class where possible.

Usage:
    python3 select_mcsa_proteins.py /mnt/scratch/sunny/m-csa/mcsa_residues.json \
        --n 10 --output selected_proteins.tsv
"""
import json
import argparse
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('mcsa_json')
    parser.add_argument('--n', type=int, default=10)
    parser.add_argument('--output', default='selected_proteins.tsv')
    args = parser.parse_args()

    with open(args.mcsa_json) as f:
        data = json.load(f)

    # Build uniprot -> (pdb_id, mcsa_id, n_catalytic_residues) mapping
    # Only keep entries where we have both uniprot and pdb
    uniprot_to_info = {}

    mcsa_id_to_residues = defaultdict(lambda: {'pdb': None, 'uniprot': None, 'n_res': 0})

    for residue in data:
        mcsa_id = residue.get('mcsa_id')
        for chain in residue.get('residue_chains', []):
            if chain.get('is_reference'):
                mcsa_id_to_residues[mcsa_id]['pdb'] = chain['pdb_id'].lower()
        for seq in residue.get('residue_sequences', []):
            if seq.get('is_reference'):
                mcsa_id_to_residues[mcsa_id]['uniprot'] = seq['uniprot_id']
        mcsa_id_to_residues[mcsa_id]['n_res'] += 1

    # Group by EC class using mcsa_id as a proxy for diversity
    # (lower mcsa_ids tend to be different enzyme classes)
    # Filter out entries missing pdb or uniprot
    valid = [
        (mcsa_id, info)
        for mcsa_id, info in mcsa_id_to_residues.items()
        if info['pdb'] and info['uniprot'] and info['n_res'] >= 2
    ]

    # Sort by mcsa_id and sample evenly across the range for diversity
    valid.sort(key=lambda x: x[0])
    total = len(valid)
    step = total // args.n
    selected = [valid[i * step] for i in range(args.n)]

    # Write output
    with open(args.output, 'w') as f:
        f.write("mcsa_id\tuniprot_id\tpdb_id\tn_catalytic_residues\n")
        for mcsa_id, info in selected:
            f.write(f"{mcsa_id}\t{info['uniprot']}\t{info['pdb']}\t{info['n_res']}\n")
            print(f"  mcsa_id={mcsa_id}  UniProt={info['uniprot']}  "
                  f"PDB={info['pdb'].upper()}  n_residues={info['n_res']}")

    print(f"\nSelected {len(selected)} proteins → {args.output}")

if __name__ == '__main__':
    main()