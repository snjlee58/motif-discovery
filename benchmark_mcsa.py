#!/usr/bin/env python3
"""
Benchmark conservation-based motif predictions against M-CSA ground truth.

Uses the parsed TSV from parse_catalytic_residues_merged.py (catalytic_residues_separated_aligned.tsv).

Usage:
    python3 benchmark_mcsa.py <conservation_json> <mcsa_tsv> <mapping_json> --pdb-id 1btl --top-n 5
"""

import json
import sys
from pathlib import Path
from typing import List, Set, Dict, Tuple

def extract_mcsa_residues(mcsa_tsv: str, pdb_id: str, chain: str = None,
                          references_only: bool = True) -> Set[int]:
    """
    Extract catalytic residue positions from parsed M-CSA TSV.
    
    TSV format: mcsa_id  pdb_id  chain_residues  is_reference
    chain_residues example: A7,A70,A178,A180,A8,A147
    
    Args:
        mcsa_tsv: Path to the parsed TSV file
        pdb_id: PDB ID to search for
        chain: Chain letter to filter (default: auto-detect first chain)
        references_only: Only use is_reference=True entries
    
    Returns:
        Set of residue positions (auth_resid as integers)
    """
    import csv
    
    pdb_id_lower = pdb_id.lower()
    catalytic_positions = set()
    found_pdbs = set()
    matched_rows = []
    
    with open(mcsa_tsv, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            row_pdb = row['pdb_id'].strip().lower()
            found_pdbs.add(row_pdb)
            
            if row_pdb != pdb_id_lower:
                continue
            
            is_ref = row['is_reference'].strip().upper() == 'TRUE'
            if references_only and not is_ref:
                continue
            
            matched_rows.append(row)
    
    if not matched_rows:
        print(f"  WARNING: PDB {pdb_id_lower} not found in M-CSA TSV"
              f" (references_only={references_only})")
        print(f"  Available PDBs (first 20): {sorted(found_pdbs)[:20]}")
        return catalytic_positions
    
    print(f"  Found {len(matched_rows)} matching row(s) for PDB {pdb_id_lower}")
    
    # Parse chain_residues from matched rows
    for row in matched_rows:
        residues_str = row['residues'].strip()
        for token in residues_str.split(','):
            token = token.strip()
            if token == '-' or not token:
                continue
            
            # Parse chain letter + residue number (e.g. "A70")
            res_chain = token[0]
            try:
                res_num = int(token[1:])
            except ValueError:
                print(f"  WARNING: Could not parse residue token: {token}")
                continue
            
            # Filter by chain if specified
            if chain and res_chain != chain.upper():
                continue
            
            catalytic_positions.add(res_num)
    
    # If no chain specified but multiple chains found, warn about it
    if not chain and len(matched_rows) > 1:
        print(f"  NOTE: Multiple chains found. Using residues from all chains.")
        print(f"  Specify --chain to restrict to one chain.")
    
    return catalytic_positions

def get_top_conserved_positions(conservation_data: Dict, alignment_mapping: Dict,
                                top_n: int = 5,
                                exclude_gaps: bool = True,
                                min_identity: float = 0.0) -> List[int]:
    """
    Extract top N conserved positions mapped to PDB residue IDs.
    
    Args:
        conservation_data: Conservation JSON
        alignment_mapping: Dict mapping alignment_column -> resid
        top_n: Number of top positions
        exclude_gaps: Filter high-gap positions
        min_identity: Minimum identity threshold
    
    Returns:
        List of PDB residue IDs (resid)
    """
    positions = conservation_data['positions']
    
    # Filter
    filtered = []
    for pos in positions:
        if exclude_gaps and pos['gap_frequency'] > 0.5:
            continue
        if pos['identity'] < min_identity:
            continue
        if pos['consensus'] in ['X', '-']:
            continue
        
        # Map alignment column to resid
        aln_col = pos['position']
        resid = alignment_mapping.get(str(aln_col))  # JSON keys are strings
        
        if resid is None:
            # This alignment column has a gap in the query sequence
            continue
        
        pos['resid'] = resid  # Add for later reference
        filtered.append(pos)
    
    # Sort by conservation score
    sorted_positions = sorted(filtered, key=lambda p: p['conservation'], reverse=True)
    
    # Return resid values
    return [p['resid'] for p in sorted_positions[:top_n]]

def calculate_metrics(predicted: List[int], true_set: Set[int]) -> Dict[str, float]:
    """
    Calculate precision, recall, and F1 score.
    
    Args:
        predicted: List of predicted positions
        true_set: Set of true catalytic positions
    
    Returns:
        Dict with precision, recall, f1, tp, fp, fn counts
    """
    predicted_set = set(predicted)
    
    tp = len(predicted_set & true_set)  # True positives
    fp = len(predicted_set - true_set)  # False positives
    fn = len(true_set - predicted_set)  # False negatives
    
    precision = tp / len(predicted_set) if predicted_set else 0.0
    recall = tp / len(true_set) if true_set else 0.0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0
    
    return {
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'tp': tp,
        'fp': fp,
        'fn': fn,
        'n_predicted': len(predicted_set),
        'n_true': len(true_set)
    }

def print_results(predicted: List[int], true_set: Set[int], metrics: Dict, 
                 conservation_data: Dict):
    """Print detailed benchmarking results."""
    predicted_set = set(predicted)
    
    print("\n" + "="*70)
    print("BASELINE PERFORMANCE: Top-N Conservation")
    print("="*70)
    
    print(f"\nM-CSA Ground Truth: {sorted(true_set)}")
    print(f"Predicted (Top {len(predicted)}): {predicted}")
    
    print("\n--- Detailed Analysis ---")
    
    # True Positives
    tp_positions = predicted_set & true_set
    if tp_positions:
        print(f"\n✓ TRUE POSITIVES ({len(tp_positions)}):")
        for pos in sorted(tp_positions):
            pos_data = next((p for p in conservation_data['positions'] if p['position'] == pos), None)
            if pos_data:
                print(f"   Position {pos} ({pos_data['consensus']}): "
                      f"conservation={pos_data['conservation']:.4f}")
    
    # False Positives
    fp_positions = predicted_set - true_set
    if fp_positions:
        print(f"\n✗ FALSE POSITIVES ({len(fp_positions)}):")
        for pos in sorted(fp_positions):
            pos_data = next((p for p in conservation_data['positions'] if p['position'] == pos), None)
            if pos_data:
                print(f"   Position {pos} ({pos_data['consensus']}): "
                      f"conservation={pos_data['conservation']:.4f}")
    
    # False Negatives
    fn_positions = true_set - predicted_set
    if fn_positions:
        print(f"\n✗ FALSE NEGATIVES (Missed catalytic residues: {len(fn_positions)}):")
        for pos in sorted(fn_positions):
            pos_data = next((p for p in conservation_data['positions'] if p['position'] == pos), None)
            if pos_data:
                print(f"   Position {pos} ({pos_data['consensus']}): "
                      f"conservation={pos_data['conservation']:.4f}, "
                      f"rank={get_conservation_rank(conservation_data, pos)}")
    
    print("\n--- Performance Metrics ---")
    print(f"Precision: {metrics['precision']:.3f} ({metrics['tp']}/{metrics['n_predicted']})")
    print(f"Recall:    {metrics['recall']:.3f} ({metrics['tp']}/{metrics['n_true']})")
    print(f"F1 Score:  {metrics['f1']:.3f}")
    
    print("\n" + "="*70)
    print("This is your BASELINE. More complex models must beat this.")
    print("="*70)

def get_conservation_rank(conservation_data: Dict, position: int) -> int:
    """Get the rank of a position by conservation score."""
    positions = conservation_data['positions']
    sorted_pos = sorted(positions, key=lambda p: p['conservation'], reverse=True)
    
    for rank, pos in enumerate(sorted_pos, 1):
        if pos['position'] == position:
            return rank
    return -1

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Benchmark against M-CSA')
    parser.add_argument('conservation_json', help='Conservation scores JSON')
    parser.add_argument('mcsa_tsv', help='Parsed M-CSA residues TSV (from parse_catalytic_residues_merged.py)')
    parser.add_argument('mapping_json', help='Alignment->resid mapping JSON')
    parser.add_argument('--pdb-id', required=True, help='PDB ID (e.g., 1btl)')
    parser.add_argument('--chain', default=None, help='Chain letter to filter (default: all chains)')
    parser.add_argument('--include-homologues', action='store_true', default=False,
                        help='Also include homologue entries (default: references only)')
    parser.add_argument('--top-n', type=int, default=5, help='Top N predictions')
    parser.add_argument('--exclude-gaps', action='store_true', default=True)
    parser.add_argument('--min-identity', type=float, default=0.0)
    parser.add_argument('--output', '-o', help='Save results to JSON file')
    
    args = parser.parse_args()
    
    # Load data
    print(f"Loading conservation data: {args.conservation_json}")
    with open(args.conservation_json) as f:
        conservation_data = json.load(f)
    
    print(f"Loading M-CSA data: {args.mcsa_tsv}")
    
    print(f"Loading alignment mapping: {args.mapping_json}")
    with open(args.mapping_json) as f:
        mapping_data = json.load(f)
        alignment_mapping = mapping_data['mapping']
    
    # Extract ground truth from TSV
    references_only = not args.include_homologues
    true_positions = extract_mcsa_residues(
        args.mcsa_tsv, args.pdb_id,
        chain=args.chain,
        references_only=references_only
    )
    
    if not true_positions:
        print(f"ERROR: No M-CSA residues found for PDB {args.pdb_id}")
        sys.exit(1)
    
    # Get predictions
    predicted_positions = get_top_conserved_positions(
        conservation_data,
        alignment_mapping,
        top_n=args.top_n,
        exclude_gaps=args.exclude_gaps,
        min_identity=args.min_identity
    )
    
    # Calculate metrics
    metrics = calculate_metrics(predicted_positions, true_positions)
    
    # Print results
    print_results(predicted_positions, true_positions, metrics, conservation_data)
    
    # Save to file if requested
    if args.output:
        results = {
            'pdb_id': args.pdb_id,
            'top_n': args.top_n,
            'mcsa_ground_truth': sorted(true_positions),
            'predicted': predicted_positions,
            'metrics': metrics,
            'parameters': {
                'exclude_gaps': args.exclude_gaps,
                'min_identity': args.min_identity
            }
        }
        with open(args.output, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\nSaved results to: {args.output}")

if __name__ == '__main__':
    main()