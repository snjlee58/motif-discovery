#!/usr/bin/env python3
"""
Extract top N most conserved positions from conservation scoring results.
This is your baseline for M-CSA benchmarking.

Usage:
    python3 src/extract_top_conserved.py <conservation_json> --top-n 5
"""

import json
import sys
from pathlib import Path
from typing import List, Dict

def load_conservation_data(json_path: str) -> Dict:
    """Load conservation scores from JSON file."""
    with open(json_path, 'r') as f:
        return json.load(f)

def get_top_conserved_positions(data: Dict, top_n: int = 5, 
                                exclude_gaps: bool = True,
                                min_identity: float = 0.0,
                                alignment_mapping: Dict = None) -> List[Dict]:
    """
    Extract top N most conserved positions.
    
    Args:
        data: Conservation data from JSON
        top_n: Number of top positions to return
        exclude_gaps: Skip positions with high gap frequency (>50%)
        min_identity: Minimum identity threshold (0-1)
        alignment_mapping: Optional dict mapping alignment_column -> resid
    
    Returns:
        List of top conserved positions sorted by conservation score
    """
    positions = data['positions']
    
    # Filter if needed
    filtered = []
    for pos in positions:
        # Skip if too many gaps
        if exclude_gaps and pos['gap_frequency'] > 0.5:
            continue
        
        # Skip if below identity threshold
        if pos['identity'] < min_identity:
            continue
            
        # Skip if consensus is gap
        if pos['consensus'] == 'X' or pos['consensus'] == '-':
            continue
        
        # Add resid if mapping provided
        if alignment_mapping:
            aln_col = pos['position']
            resid = alignment_mapping.get(str(aln_col))
            if resid is None:
                # Gap in query sequence at this position
                continue
            pos['resid'] = resid
        
        filtered.append(pos)
    
    # Sort by conservation score (highest first)
    sorted_positions = sorted(filtered, key=lambda p: p['conservation'], reverse=True)
    
    return sorted_positions[:top_n]

def print_top_positions(positions: List[Dict]):
    """Print top positions in a readable format."""
    has_resid = 'resid' in positions[0] if positions else False
    
    if has_resid:
        print(f"\n{'Rank':<6} {'Col':<6} {'ResID':<7} {'AA':<4} {'Conservation':<14} {'Entropy':<12} {'Identity':<10} {'Gap%':<8}")
        print("=" * 80)
        
        for i, pos in enumerate(positions, 1):
            print(f"{i:<6} {pos['position']:<6} {pos['resid']:<7} {pos['consensus']:<4} "
                  f"{pos['conservation']:<14.6f} {pos['entropy']:<12.6f} "
                  f"{pos['identity']:<10.6f} {pos['gap_frequency']*100:<7.2f}%")
    else:
        print(f"\n{'Rank':<6} {'Pos':<6} {'AA':<4} {'Conservation':<14} {'Entropy':<12} {'Identity':<10} {'Gap%':<8}")
        print("=" * 70)
        
        for i, pos in enumerate(positions, 1):
            print(f"{i:<6} {pos['position']:<6} {pos['consensus']:<4} "
                  f"{pos['conservation']:<14.6f} {pos['entropy']:<12.6f} "
                  f"{pos['identity']:<10.6f} {pos['gap_frequency']*100:<7.2f}%")

def save_motif_residues(positions: List[Dict], output_path: str):
    """Save motif as residue list (for FoldDisco later)."""
    has_resid = 'resid' in positions[0] if positions else False
    
    with open(output_path, 'w') as f:
        if has_resid:
            residue_nums = [str(p['resid']) for p in positions]
        else:
            residue_nums = [str(p['position']) for p in positions]
        f.write(','.join(residue_nums) + '\n')
    
    id_type = "residue IDs" if has_resid else "alignment columns"
    print(f"\nSaved motif {id_type} to: {output_path}")

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Extract top conserved positions')
    parser.add_argument('conservation_json', help='Path to conservation JSON file')
    parser.add_argument('--top-n', type=int, default=5, 
                       help='Number of top positions to extract (default: 5)')
    parser.add_argument('--exclude-gaps', action='store_true', default=True,
                       help='Exclude positions with >50%% gaps')
    parser.add_argument('--min-identity', type=float, default=0.0,
                       help='Minimum identity threshold (0-1, default: 0)')
    parser.add_argument('--mapping', help='Alignment->resid mapping JSON (optional)')
    parser.add_argument('--output', '-o', help='Save residue list to file (for FoldDisco)')
    
    args = parser.parse_args()
    
    # Load data
    print(f"Loading conservation data from: {args.conservation_json}")
    data = load_conservation_data(args.conservation_json)
    
    # Load mapping if provided
    alignment_mapping = None
    if args.mapping:
        print(f"Loading alignment mapping from: {args.mapping}")
        with open(args.mapping) as f:
            mapping_data = json.load(f)
            alignment_mapping = mapping_data['mapping']
    
    print(f"\nAlignment info:")
    print(f"  Length: {data['alignment_length']} columns")
    print(f"  Sequences: {data['n_sequences']}")
    print(f"  Mean conservation: {data.get('mean_conservation', 'N/A')}")
    
    # Extract top positions
    top_positions = get_top_conserved_positions(
        data, 
        top_n=args.top_n,
        exclude_gaps=args.exclude_gaps,
        min_identity=args.min_identity,
        alignment_mapping=alignment_mapping
    )
    
    print(f"\nTop {len(top_positions)} most conserved positions:")
    print_top_positions(top_positions)
    
    # Save if requested
    if args.output:
        save_motif_residues(top_positions, args.output)
    
    # Print summary for M-CSA benchmarking
    has_resid = 'resid' in top_positions[0] if top_positions else False
    
    print("\n" + "="*70)
    print("FOR M-CSA BENCHMARKING:")
    print("="*70)
    
    if has_resid:
        residue_list = [f"{p['consensus']}{p['resid']}" for p in top_positions]
        print(f"Predicted catalytic residues: {', '.join(residue_list)}")
        print(f"Residue IDs (resid): {[p['resid'] for p in top_positions]}")
        print(f"Alignment columns: {[p['position'] for p in top_positions]}")
    else:
        residue_list = [f"{p['consensus']}{p['position']}" for p in top_positions]
        print(f"Predicted catalytic residues: {', '.join(residue_list)}")
        print(f"Alignment columns: {[p['position'] for p in top_positions]}")
        print("\nNote: Use --mapping to convert to PDB residue IDs for M-CSA comparison")
    
    print("\nNext step: Compare these to M-CSA annotations for this enzyme family")

if __name__ == '__main__':
    main()