#!/usr/bin/env python3

"""
Parse catalytic_residues_homologues.json into TSV format 
Output format: mcsa_id\tpdb_id\tchain_residues\tis_reference
Keep separate lines for each MCSA-PDB combination with aligned residue strings
"""

import json
import sys
from collections import defaultdict, OrderedDict

def parse_json_to_separate_tsv(json_file, output_file):
    """Parse the JSON file keeping separate lines for each MCSA ID-PDB ID combination"""
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    print(f"Loaded {len(data)} MCSA entries from {json_file}")
    
    # Dictionary to collect data per MCSA-PDB combination
    mcsa_pdb_data = defaultdict(lambda: {
        'chain_residues': [],
        'is_reference': False
    })
    
    # Process all entries
    for entry in data:
        mcsa_id = entry.get('mcsa_id', 'unknown')
        
        # Process residue_chains section
        if 'residue_chains' in entry:
            for residue in entry['residue_chains']:
                pdb_id = residue.get('pdb_id', '').lower()
                chain_name = residue.get('chain_name', '')
                resid = residue.get('auth_resid') or residue.get('resid', '')
                is_reference = residue.get('is_reference', False)
                
                # Skip entries without valid PDB ID
                if not pdb_id or pdb_id == '':
                    continue
                
                # Create key for MCSA-PDB combination
                key = (str(mcsa_id), pdb_id)
                
                # Set reference status if any residue is reference
                if is_reference:
                    mcsa_pdb_data[key]['is_reference'] = True
                
                # Create chain_residue string
                if resid is not None and resid != '' and chain_name:
                    chain_residue = f"{chain_name}{resid}"
                    mcsa_pdb_data[key]['chain_residues'].append(chain_residue)
                else:
                    mcsa_pdb_data[key]['chain_residues'].append('-')
    
    # Find the maximum number of residues for alignment
    max_residues = 0
    for key, info in mcsa_pdb_data.items():
        if len(info['chain_residues']) > max_residues:
            max_residues = len(info['chain_residues'])
    
    print(f"Maximum residues per entry: {max_residues}")
    
    # Write output file
    with open(output_file, 'w') as f:
        # Write header
        f.write("mcsa_id\tpdb_id\tchain_residues\tis_reference\n")
        
        total_entries = 0
        for (mcsa_id, pdb_id), info in sorted(mcsa_pdb_data.items()):
            # Remove duplicates while preserving order
            chain_residues_unique = list(dict.fromkeys(info['chain_residues']))
            
            # Pad with gaps to match maximum length
            while len(chain_residues_unique) < max_residues:
                chain_residues_unique.append('-')
            
            # Join chain residues
            chain_residues_str = ','.join(chain_residues_unique)
            
            # Write the line
            f.write(f"{mcsa_id}\t{pdb_id}\t{chain_residues_str}\t{info['is_reference']}\n")
            total_entries += 1
        
        print(f"Generated {total_entries} separate MCSA-PDB entries in {output_file}")

def parse_json_to_mcsa_merged_tsv(json_file, output_file):
    """Parse the JSON file and convert to MCSA-based merged TSV format (one line per MCSA ID)"""
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    print(f"Loaded {len(data)} MCSA entries from {json_file}")
    
    with open(output_file, 'w') as f:
        # Write header
        f.write("mcsa_id\tchain_residues\tis_reference\n")
        
        for entry in data:
            mcsa_id = entry.get('mcsa_id', 'unknown')
            
            # Collect all residues for this MCSA ID
            all_residues = []
            reference_status = []
            
            # Process residue_chains section which contains PDB information
            if 'residue_chains' in entry:
                for residue in entry['residue_chains']:
                    pdb_id = residue.get('pdb_id', '').lower()
                    chain_name = residue.get('chain_name', '')
                    resid = residue.get('auth_resid') or residue.get('resid', '')
                    is_reference = residue.get('is_reference', False)
                    
                    # Handle missing values
                    if not pdb_id or pdb_id == '':
                        pdb_id = '-'
                    if resid is None or resid == '':
                        resid = '-'
                    if not chain_name:
                        chain_name = '-'
                    
                    # Create chain_residue string
                    if resid != '-' and chain_name != '-':
                        chain_residue = f"{pdb_id}:{chain_name}{resid}"
                    else:
                        chain_residue = f"{pdb_id}:-"
                    
                    all_residues.append(chain_residue)
                    reference_status.append(str(is_reference))
            
            # If no residue_chains, check residue_sequences (fallback)
            if not all_residues and 'residue_sequences' in entry:
                for residue in entry['residue_sequences']:
                    uniprot_id = residue.get('uniprot_id', '-')
                    resid = residue.get('resid', '')
                    is_reference = residue.get('is_reference', False)
                    
                    if resid is None or resid == '':
                        resid = '-'
                    
                    chain_residue = f"{uniprot_id}:{resid}"
                    all_residues.append(chain_residue)
                    reference_status.append(str(is_reference))
            
            # If still no residues, add a placeholder
            if not all_residues:
                all_residues = ['-']
                reference_status = ['False']
            
            # Join all residues and reference statuses
            chain_residues_merged = ','.join(all_residues)
            is_reference_merged = ','.join(reference_status)
            
            # Write the merged line
            f.write(f"{mcsa_id}\t{chain_residues_merged}\t{is_reference_merged}\n")
        
        print(f"Generated {len(data)} MCSA-merged entries in {output_file}")

def main():
    json_file = "catalytic_residues_homologues.json"
    
    # Generate MCSA-based merged version (original functionality)
    output_file1 = "catalytic_residues_merged_mcsa.tsv"
    parse_json_to_mcsa_merged_tsv(json_file, output_file1)
    
    # Generate separate MCSA-PDB version with aligned residues
    output_file2 = "catalytic_residues_separated_aligned.tsv"
    parse_json_to_separate_tsv(json_file, output_file2)
    
    print(f"\nParsing complete!")
    print(f"Input: {json_file}")
    print(f"MCSA-merged output: {output_file1}")
    print(f"Separated-aligned output: {output_file2}")
    
    # Show first few lines of both outputs
    for output_file in [output_file1, output_file2]:
        print(f"\nFirst 5 lines of {output_file}:")
        try:
            with open(output_file, 'r') as f:
                for i, line in enumerate(f):
                    if i >= 5:
                        break
                    print(line.strip())
        except FileNotFoundError:
            print(f"File {output_file} not found")

if __name__ == "__main__":
    main()