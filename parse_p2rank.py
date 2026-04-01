#!/usr/bin/env python3
"""
Parse P2Rank residue-level output into JSON format for pipeline integration.

P2Rank outputs a *_residues.csv with per-residue binding site scores.
This script extracts those scores keyed by auth_resid.

Usage:
    python3 parse_p2rank.py <p2rank_residues.csv> -o p2rank_scores.json [--chain A]
"""

import csv
import json
import argparse
import re


def parse_p2rank_residues(residues_csv: str, chain: str = None):
    """
    Parse P2Rank residues CSV file.
    
    P2Rank residues.csv columns:
      chain, residue_label, residue_name, score, zscore, probability, 
      pocket, sas, buried
    
    residue_label format: "A_70" (chain_resnum) or just "70"
    
    Args:
        residues_csv: Path to P2Rank *_residues.csv
        chain: Chain to filter (default: first chain or all)
    
    Returns:
        Dict of {auth_resid: {"score": float, "probability": float, "pocket": int}}
    """
    scores = {}
    target_chain = chain.upper() if chain else None
    
    with open(residues_csv, 'r') as f:
        reader = csv.DictReader(f)
        
        # Clean up column names (P2Rank adds spaces)
        reader.fieldnames = [name.strip() for name in reader.fieldnames]
        
        for row in reader:
            # Clean whitespace from values
            row = {k: v.strip() if isinstance(v, str) else v for k, v in row.items()}
            
            # Parse chain and residue number
            res_chain = row.get('chain', '').strip()
            res_label = row.get('residue_label', '').strip()
            
            # Auto-detect first chain
            if target_chain is None and res_chain:
                target_chain = res_chain
            
            # Filter by chain
            if target_chain and res_chain and res_chain != target_chain:
                continue
            
            # Extract residue number from label
            # Format can be "A_70", "70", or "A 70"
            resnum = None
            if '_' in res_label:
                parts = res_label.split('_')
                try:
                    resnum = int(parts[-1])
                except ValueError:
                    pass
            else:
                try:
                    resnum = int(res_label)
                except ValueError:
                    # Try extracting digits
                    digits = re.findall(r'\d+', res_label)
                    if digits:
                        resnum = int(digits[-1])
            
            if resnum is None:
                continue
            
            # Extract scores
            try:
                score = float(row.get('score', 0))
                probability = float(row.get('probability', 0))
                pocket = int(row.get('pocket', 0)) if row.get('pocket', '0') else 0
            except (ValueError, TypeError):
                score = 0.0
                probability = 0.0
                pocket = 0
            
            scores[resnum] = {
                'score': score,
                'probability': probability,
                'pocket': pocket,
                'residue_name': row.get('residue_name', '').strip(),
            }
    
    return scores


def main():
    parser = argparse.ArgumentParser(description='Parse P2Rank output to JSON')
    parser.add_argument('residues_csv', help='P2Rank *_residues.csv file')
    parser.add_argument('-o', '--output', default='p2rank_scores.json', help='Output JSON')
    parser.add_argument('--chain', default=None, help='Chain to extract (default: first)')
    
    args = parser.parse_args()
    
    print(f"Parsing P2Rank output: {args.residues_csv}")
    scores = parse_p2rank_residues(args.residues_csv, chain=args.chain)
    
    # Save
    output_data = {
        'source': 'p2rank',
        'input': args.residues_csv,
        'chain': args.chain or 'auto',
        'n_residues': len(scores),
        'residues': {str(k): v for k, v in scores.items()}
    }
    
    with open(args.output, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"Extracted scores for {len(scores)} residues")
    print(f"Saved to: {args.output}")
    
    # Summary
    if scores:
        probs = [v['probability'] for v in scores.values()]
        in_pocket = sum(1 for v in scores.values() if v['pocket'] > 0)
        high_prob = sum(1 for p in probs if p > 0.5)
        print(f"Residues in a predicted pocket: {in_pocket}")
        print(f"Residues with probability > 0.5: {high_prob}")
        print(f"Probability range: {min(probs):.3f} - {max(probs):.3f}")


if __name__ == '__main__':
    main()
