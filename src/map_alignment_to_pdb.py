#!/usr/bin/env python3
"""
Map alignment column numbers to PDB residue numbers.

This creates the mapping needed for benchmarking:
  alignment_column -> PDB auth_resid -> compare with M-CSA

Usage:
    python3 src/map_alignment_to_pdb.py <foldmason_msa.fa> <pdb_id> -o mapping.json
    python3 src/map_alignment_to_pdb.py <foldmason_msa.fa> <pdb_id> --pdb-file 1BTL.pdb -o mapping.json
"""

import json
from Bio import AlignIO
from pathlib import Path


def extract_auth_resids_from_pdb(pdb_file: str, chain: str = None):
    """
    Extract ordered list of (chain, auth_resid) from a PDB file.
    Only considers ATOM records (not HETATM). Returns unique residues 
    in the order they appear, handling insertion codes.
    
    Args:
        pdb_file: Path to PDB file
        chain: Chain letter to filter (default: first chain found)
    
    Returns:
        List of auth_resid integers in residue order
    """
    seen = set()
    residues = []
    target_chain = chain.upper() if chain else None
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            
            ch = line[21].strip()
            
            # Auto-detect first chain if not specified
            if target_chain is None:
                target_chain = ch
            
            if ch != target_chain:
                continue
            
            # PDB format: columns 22-26 = resSeq, column 26 = insertion code
            resseq_str = line[22:26].strip()
            icode = line[26].strip()
            
            try:
                resseq = int(resseq_str)
            except ValueError:
                continue
            
            # Use (resseq, icode) as unique key to handle insertion codes
            key = (resseq, icode)
            if key not in seen:
                seen.add(key)
                residues.append(resseq)
    
    print(f"  Extracted {len(residues)} residues from PDB file (chain {target_chain})")
    if residues:
        print(f"  Auth resid range: {residues[0]} - {residues[-1]}")
    
    return residues

def create_alignment_to_resid_mapping(alignment_file: str, pdb_id: str, 
                                      pdb_file: str = None, chain: str = None,
                                      uniprot_id: str = None):
    """
    Create mapping from alignment columns to PDB residue IDs.
    
    If pdb_file is provided, uses auth_resid from the PDB file.
    Otherwise falls back to sequential numbering (1, 2, 3...).
    
    Args:
        alignment_file: Path to FoldMason MSA file
        pdb_id: PDB ID to find in alignment (e.g., '1BTL' or '1btl')
        pdb_file: Optional path to PDB file for auth_resid extraction
        chain: Chain letter for PDB parsing (default: first chain)
    
    Returns:
        dict: {alignment_column: resid}
    """
    # Load alignment
    alignment = AlignIO.read(alignment_file, "fasta")
    
    # Find the query sequence - try PDB ID first, then UniProt/AF ID as fallback
    pdb_id_upper = pdb_id.upper()
    query_seq = None
    query_id = None
    
    for record in alignment:
        if pdb_id_upper in record.id.upper():
            query_seq = str(record.seq)
            query_id = record.id
            break
    
    # Fallback: search for AlphaFold entry matching UniProt ID
    if query_seq is None and uniprot_id:
        uniprot_upper = uniprot_id.upper()
        for record in alignment:
            if uniprot_upper in record.id.upper():
                query_seq = str(record.seq)
                query_id = record.id
                print(f"PDB {pdb_id} not found in alignment, using AlphaFold entry: {query_id}")
                break
    
    if query_seq is None:
        available = [r.id for r in alignment[:10]]
        raise ValueError(f"Could not find {pdb_id} in alignment. Available IDs: {available}")
    
    print(f"Found query sequence: {query_id}")
    print(f"Alignment length: {len(query_seq)} columns")
    
    # Get auth_resid list from PDB file, or fall back to sequential
    auth_resids = None
    if pdb_file and Path(pdb_file).exists():
        print(f"Reading auth_resid from PDB file: {pdb_file}")
        auth_resids = extract_auth_resids_from_pdb(pdb_file, chain=chain)
        
        # Count non-gap residues in alignment
        n_residues_in_aln = sum(1 for aa in query_seq if aa != '-')
        if len(auth_resids) != n_residues_in_aln:
            print(f"  WARNING: PDB has {len(auth_resids)} residues but alignment "
                  f"has {n_residues_in_aln} non-gap positions")
            if abs(len(auth_resids) - n_residues_in_aln) > 5:
                print(f"  WARNING: Large mismatch! Falling back to sequential numbering.")
                auth_resids = None
    else:
        if pdb_file:
            print(f"WARNING: PDB file not found at {pdb_file}, using sequential numbering")
        else:
            print("No PDB file provided, using sequential residue numbering")
    
    # Create mapping: alignment_column -> resid
    mapping = {}
    residue_idx = 0  # Index into auth_resids list
    
    for col_idx, aa in enumerate(query_seq, start=1):
        if aa != '-':  # Not a gap
            if auth_resids and residue_idx < len(auth_resids):
                mapping[col_idx] = auth_resids[residue_idx]
            else:
                mapping[col_idx] = residue_idx + 1  # Sequential fallback
            residue_idx += 1
    
    n_residues = residue_idx
    print(f"Total residues in {pdb_id}: {n_residues}")
    print(f"Mapped {len(mapping)} alignment columns to residue IDs")
    if auth_resids:
        print(f"Using auth_resid numbering from PDB file")
    else:
        print(f"Using sequential numbering (1-based)")
    
    return mapping, query_id, n_residues

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Map alignment columns to residue IDs')
    parser.add_argument('alignment', help='FoldMason MSA file (foldmason_result_aa.fa)')
    parser.add_argument('pdb_id', help='PDB ID to extract (e.g., 1BTL)')
    parser.add_argument('--pdb-file', help='PDB file for auth_resid extraction (recommended)')
    parser.add_argument('--uniprot', default=None, help='UniProt ID as fallback for finding query in alignment')
    parser.add_argument('--chain', default=None, help='Chain letter in PDB file (default: first chain)')
    parser.add_argument('-o', '--output', help='Output JSON file', default='alignment_mapping.json')
    
    args = parser.parse_args()
    
    # Create mapping
    mapping, query_id, n_residues = create_alignment_to_resid_mapping(
        args.alignment, 
        args.pdb_id,
        pdb_file=args.pdb_file,
        chain=args.chain,
        uniprot_id=args.uniprot
    )
    
    # Save to file
    output_data = {
        'pdb_id': args.pdb_id,
        'query_id': query_id,
        'n_residues': n_residues,
        'alignment_length': max(mapping.keys()),
        'mapping': {str(k): v for k, v in mapping.items()}  # JSON keys must be strings
    }
    
    with open(args.output, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\nSaved mapping to: {args.output}")
    
    # Show first few mappings
    print("\nFirst 10 mappings (alignment_column -> resid):")
    for col in sorted(mapping.keys())[:10]:
        print(f"  Column {col} -> Residue {mapping[col]}")

if __name__ == '__main__':
    main()