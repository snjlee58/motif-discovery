#!/usr/bin/env python3
"""
Parse M-CSA representative (non-homologue) dataset from the API into TSV format.

Input:  catalytic_residues.json (from https://www.ebi.ac.uk/thornton-srv/m-csa/api/residues/?format=json)
Output: catalytic_residues_parsed.tsv with columns: mcsa_id, pdb_id, residues, is_reference

Each row is one MCSA-PDB combination with all catalytic residues comma-separated
as ChainResidue (e.g., A7,A70,A178).

Usage:
    python3 parse_mcsa_representative.py catalytic_residues_representative_raw.json -o catalytic_residues_representative_parsed.tsv
"""

import json
import argparse
from collections import defaultdict


def parse_mcsa_representative(json_file, output_file):
    with open(json_file, 'r') as f:
        data = json.load(f)

    print(f"Loaded {len(data)} residue entries from {json_file}")

    # Group by (mcsa_id, pdb_id)
    # Each JSON entry is one residue, so we need to collect all residues
    # for the same mcsa_id + pdb_id combination
    groups = defaultdict(lambda: {
        'residues': [],
        'is_reference': False
    })

    skipped = 0
    for entry in data:
        mcsa_id = entry.get('mcsa_id', 'unknown')

        # Each entry has residue_chains with PDB-level info
        if 'residue_chains' not in entry:
            skipped += 1
            continue

        for res in entry['residue_chains']:
            pdb_id = res.get('pdb_id', '').lower()
            chain_name = res.get('chain_name', '')
            # Prefer auth_resid (matches PDB file numbering), fall back to resid
            auth_resid = res.get('auth_resid')
            resid = auth_resid if auth_resid is not None else res.get('resid', '')
            is_reference = res.get('is_reference', False)

            if not pdb_id or not chain_name or resid == '' or resid is None:
                skipped += 1
                continue

            key = (str(mcsa_id), pdb_id)

            if is_reference:
                groups[key]['is_reference'] = True

            chain_residue = f"{chain_name}{resid}"
            # Avoid duplicates (same residue can appear in multiple roles)
            if chain_residue not in groups[key]['residues']:
                groups[key]['residues'].append(chain_residue)

    print(f"Skipped {skipped} entries with missing data")
    print(f"Found {len(groups)} unique MCSA-PDB combinations")

    # Count references vs homologues
    n_ref = sum(1 for g in groups.values() if g['is_reference'])
    n_hom = len(groups) - n_ref
    print(f"  Reference entries: {n_ref}")
    print(f"  Non-reference entries: {n_hom}")

    # Stats on catalytic residues per entry
    n_res_list = [len(g['residues']) for g in groups.values()]
    if n_res_list:
        print(f"  Catalytic residues per entry: min={min(n_res_list)}, "
              f"max={max(n_res_list)}, mean={sum(n_res_list)/len(n_res_list):.1f}")

    # Unique PDB IDs
    unique_pdbs = set(pdb_id for (_, pdb_id) in groups.keys())
    print(f"  Unique PDB IDs: {len(unique_pdbs)}")

    # Unique MCSA IDs
    unique_mcsa = set(mcsa_id for (mcsa_id, _) in groups.keys())
    print(f"  Unique MCSA IDs: {len(unique_mcsa)}")

    # Write output
    with open(output_file, 'w') as f:
        f.write("mcsa_id\tpdb_id\tresidues\tis_reference\n")

        for (mcsa_id, pdb_id), info in sorted(groups.items(),
                key=lambda x: (int(x[0][0]) if x[0][0].isdigit() else 0, x[0][1])):
            residues_str = ','.join(info['residues'])
            f.write(f"{mcsa_id}\t{pdb_id}\t{residues_str}\t{info['is_reference']}\n")

    print(f"\nWritten to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Parse M-CSA representative dataset JSON into TSV'
    )
    parser.add_argument('json_file', help='Path to M-CSA residues JSON file')
    parser.add_argument('-o', '--output', default='catalytic_residues_representative_parsed.tsv',
                        help='Output TSV file')
    args = parser.parse_args()

    parse_mcsa_representative(args.json_file, args.output)


if __name__ == '__main__':
    main()
