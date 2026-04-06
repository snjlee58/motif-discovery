#!/usr/bin/env python3
"""
Benchmark conservation-based motif predictions against M-CSA ground truth.

Uses the parsed TSV from parse_catalytic_residues_merged.py (catalytic_residues_homologues_parsed.tsv).

Usage:
    python3 benchmark_mcsa.py <conservation_json> <mcsa_tsv> <mapping_json> --pdb-id 1btl --top-n 5
"""

import json
import math
import sys
from pathlib import Path
from typing import List, Set, Dict, Tuple


def parse_ca_coordinates(pdb_file: str) -> Dict[int, Tuple[float, float, float]]:
    """
    Parse CA (alpha carbon) coordinates from a PDB file.
    
    Returns:
        Dict mapping auth_resid -> (x, y, z) coordinates
    """
    coords = {}
    with open(pdb_file) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            try:
                resid = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                # Keep first occurrence (in case of alt conformations)
                if resid not in coords:
                    coords[resid] = (x, y, z)
            except (ValueError, IndexError):
                continue
    return coords


def compute_spatial_clustering(candidates: List[Dict], ca_coords: Dict[int, Tuple],
                               radius: float = 15.0, top_fraction: int = 3) -> Dict[int, float]:
    """
    Compute spatial clustering score for candidate residues.
    
    For each candidate, count how many other high-scoring candidates are
    within `radius` Angstroms. Normalize to [0, 1].
    
    This is a two-pass approach:
      1. Take the top N*top_fraction candidates by preliminary score
      2. For each, count neighbors within radius among this pool
      3. Normalize: cluster_score = neighbor_count / max_neighbor_count
    
    Args:
        candidates: List of position dicts with 'resid' and 'combined_score'
        ca_coords: Dict mapping resid -> (x, y, z)
        radius: Distance cutoff in Angstroms for "nearby"
        top_fraction: Consider top N*this candidates for clustering pool
    
    Returns:
        Dict mapping resid -> clustering score [0, 1]
    """
    if not candidates or not ca_coords:
        return {}
    
    # Get coordinates for candidates that have them
    candidates_with_coords = []
    for c in candidates:
        resid = c['resid']
        if resid in ca_coords:
            candidates_with_coords.append((resid, ca_coords[resid]))
    
    if len(candidates_with_coords) < 2:
        return {}
    
    # Count neighbors within radius for each candidate
    neighbor_counts = {}
    for i, (res_i, coord_i) in enumerate(candidates_with_coords):
        count = 0
        for j, (res_j, coord_j) in enumerate(candidates_with_coords):
            if i == j:
                continue
            dist = math.sqrt(
                (coord_i[0] - coord_j[0]) ** 2 +
                (coord_i[1] - coord_j[1]) ** 2 +
                (coord_i[2] - coord_j[2]) ** 2
            )
            if dist <= radius:
                count += 1
        neighbor_counts[res_i] = count
    
    # Normalize to [0, 1]
    max_count = max(neighbor_counts.values()) if neighbor_counts else 1
    if max_count == 0:
        return {r: 0.0 for r in neighbor_counts}
    
    return {r: count / max_count for r, count in neighbor_counts.items()}


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
        # Support both column names (depends on parse script version)
        residues_str = (row.get('chain_residues') or row.get('residues') or '').strip()
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
                                conservation_threshold: float = None,
                                exclude_gaps: bool = True,
                                min_identity: float = 0.0,
                                exclude_structural: bool = False,
                                use_catalytic_propensity: bool = False,
                                p2rank_scores: Dict = None,
                                pdb_file: str = None) -> List[int]:
    """
    Extract top conserved positions mapped to PDB residue IDs.
    
    Uses either top-N or conservation threshold for selection.
    Combined score = w1*conservation + w2*p2rank + w3*propensity + w4*clustering
    
    Two-pass approach:
      Pass 1: Score with conservation + propensity + p2rank
      Pass 2: Compute spatial clustering among top candidates, add as bonus, re-rank
    
    Args:
        conservation_data: Conservation JSON
        alignment_mapping: Dict mapping alignment_column -> resid
        top_n: Number of top positions (ignored if threshold is set)
        conservation_threshold: Min conservation score (overrides top_n)
        exclude_gaps: Filter high-gap positions
        min_identity: Minimum identity threshold
        exclude_structural: Hard-filter G, P, A (aggressive, not recommended)
        use_catalytic_propensity: Weight scores by M-CSA catalytic propensity
        p2rank_scores: Dict of {auth_resid: {"score": float, "probability": float, ...}}
                       from parse_p2rank.py output
        pdb_file: Path to PDB file for spatial clustering (optional)
    
    Returns:
        List of PDB residue IDs (resid)
    """
    # Catalytic propensity from M-CSA (Ribeiro et al. 2017)
    CATALYTIC_PROPENSITY = {
        'H': 8.01, # Histidine (His)
        'C': 4.66, # Cysteine (Cys)
        'D': 3.04, # Aspartic acid (Asp)
        'E': 2.09, # Glutamic acid (Glu)
        'K': 1.58, # Lysine (Lys)
        'R': 1.81, # Arginine (Arg)
        'S': 0.95, # Serine (Ser)
        'T': 0.56, # Threonine (Thr)
        'Y': 1.84, # Tyrosine (Tyr)
        'N': 0.97, # Asparagine (Asn)
        'Q': 0.48, # Glutamine (Gln)
        'W': 1.03, # Tryptophan (Trp)
        'F': 0.40, # Phenylalanine (Phe)
        'M': 0.26, # Methionine (Met)
        'I': 0.04, # Isoleucine (Ile)
        'L': 0.04, # Leucine (Leu)
        'V': 0.03, # Valine (Val)
        'A': 0.01, # Alanine (Ala)
        'G': 0.02, # Glycine (Gly)
        'P': 0.03, # Proline (Pro)
        'X': 0.5, # Unknown
    }
    
    STRUCTURAL_AAS = {'G', 'P', 'A'}
    
    positions = conservation_data['positions']
    
    # Filter and score
    filtered = []
    n_structural_skipped = 0
    n_p2rank_matched = 0
    
    for pos in positions:
        if exclude_gaps and pos['gap_frequency'] > 0.5:
            continue
        if pos['identity'] < min_identity:
            continue
        if pos['consensus'] in ['X', '-']:
            continue
        
        if exclude_structural and pos['consensus'] in STRUCTURAL_AAS:
            n_structural_skipped += 1
            continue
        
        # Map alignment column to resid
        aln_col = pos['position']
        resid = alignment_mapping.get(str(aln_col))
        
        if resid is None:
            continue
        
        pos['resid'] = resid
        
        # === ADDITIVE SCORING ===
        # Each signal normalized to [0, 1], combined additively
        
        # Signal 1: Conservation (already 0-1)
        s_conservation = pos['conservation']
        
        # Signal 2: Catalytic propensity
        s_propensity = 0.0
        if use_catalytic_propensity:
            aa = pos['consensus'].upper()
            propensity = CATALYTIC_PROPENSITY.get(aa, 0.5)
            # Normalize: max propensity is ~8.01 (His)
            s_propensity = min(propensity / 8.01, 1.0)
            pos['catalytic_propensity'] = propensity
        
        # Signal 3: P2Rank binding site probability
        s_p2rank = 0.0
        if p2rank_scores:
            p2rank_data = p2rank_scores.get(str(resid)) or p2rank_scores.get(resid)
            if p2rank_data:
                n_p2rank_matched += 1
                if isinstance(p2rank_data, dict):
                    p2rank_prob = p2rank_data.get('probability', 0.0)
                else:
                    p2rank_prob = float(p2rank_data)
                s_p2rank = p2rank_prob  # already 0-1
                pos['p2rank_probability'] = p2rank_prob
            else:
                # No P2Rank data for this residue — leave at 0, don't penalize
                pos['p2rank_probability'] = None
        
        # Weighted additive combination
        # Conservation is the backbone, other signals boost
        W_CONSERVATION = 1.0
        W_P2RANK = 0.35
        W_PROPENSITY = 0.25
        
        score = (W_CONSERVATION * s_conservation
                 + W_P2RANK * s_p2rank
                 + W_PROPENSITY * s_propensity)
        
        pos['combined_score'] = score
        filtered.append(pos)
    
    if exclude_structural and n_structural_skipped > 0:
        print(f"  Excluded {n_structural_skipped} structural residues (G, P, A)")
    if p2rank_scores:
        print(f"  P2Rank scores matched: {n_p2rank_matched}/{len(filtered)} residues")
    
    # === PASS 1: Sort by preliminary score (conservation + propensity + p2rank) ===
    sorted_pass1 = sorted(filtered, key=lambda p: p['combined_score'], reverse=True)
    
    # === PASS 2: Spatial clustering bonus ===
    W_CLUSTERING = 0.30
    
    if pdb_file and Path(pdb_file).exists():
        ca_coords = parse_ca_coordinates(pdb_file)
        
        # Take top candidates for clustering pool (3x top_n, or top 30, whichever is larger)
        pool_size = max((top_n or 10) * 3, 30)
        pool = sorted_pass1[:pool_size]
        
        clustering_scores = compute_spatial_clustering(pool, ca_coords, radius=15.0)
        
        if clustering_scores:
            n_clustered = sum(1 for s in clustering_scores.values() if s > 0)
            print(f"  Spatial clustering: {n_clustered}/{len(clustering_scores)} residues have neighbors within 15Å")
            
            # Add clustering bonus and recompute combined score
            for pos in filtered:
                resid = pos['resid']
                s_cluster = clustering_scores.get(resid, 0.0)
                pos['cluster_score'] = s_cluster
                pos['combined_score'] = pos['combined_score'] + W_CLUSTERING * s_cluster
        else:
            print(f"  Spatial clustering: no coordinates found, skipping")
            for pos in filtered:
                pos['cluster_score'] = 0.0
    else:
        if pdb_file:
            print(f"  Spatial clustering: PDB file not found ({pdb_file}), skipping")
        for pos in filtered:
            pos['cluster_score'] = 0.0
    
    # === Final sort after clustering ===
    sorted_positions = sorted(filtered, key=lambda p: p['combined_score'], reverse=True)
    
    # Log what signals are active
    signals = []
    if use_catalytic_propensity:
        signals.append("catalytic propensity")
    if p2rank_scores:
        signals.append("P2Rank binding site")
    if pdb_file and Path(pdb_file).exists():
        signals.append("spatial clustering")
    if signals:
        print(f"  Active signals: conservation + {' + '.join(signals)}")
    
    # Select by threshold or top-N
    if conservation_threshold is not None:
        selected = [p for p in sorted_positions if p['combined_score'] >= conservation_threshold]
        print(f"  Threshold {conservation_threshold}: {len(selected)} residues above cutoff")
        return [p['resid'] for p in selected]
    else:
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

def _build_resid_lookup(conservation_data: Dict) -> Dict[int, Dict]:
    """Build a lookup from auth_resid -> position data.
    
    Uses the 'resid' field that get_top_conserved_positions adds to each position.
    Falls back to 'auth_resid' if present.
    """
    lookup = {}
    for p in conservation_data['positions']:
        resid = p.get('resid') or p.get('auth_resid')
        if resid is not None:
            lookup[resid] = p
    return lookup

def print_results(predicted: List[int], true_set: Set[int], metrics: Dict, 
                 conservation_data: Dict, alignment_mapping: Dict = None):
    """Print detailed benchmarking results."""
    predicted_set = set(predicted)
    
    # Build resid -> position data lookup
    resid_lookup = _build_resid_lookup(conservation_data)
    
    # Also build from alignment_mapping as fallback
    if alignment_mapping:
        col_to_resid = {}
        resid_to_col = {}
        for col_str, resid in alignment_mapping.items():
            col_to_resid[int(col_str)] = resid
            resid_to_col[resid] = int(col_str)
        
        # Fill in any missing entries
        for p in conservation_data['positions']:
            aln_col = p['position']
            resid = col_to_resid.get(aln_col)
            if resid is not None and resid not in resid_lookup:
                resid_lookup[resid] = p
    
    def _get_pos_data(resid):
        return resid_lookup.get(resid)
    
    print("\n" + "="*70)
    print("BASELINE PERFORMANCE: Top-N Conservation")
    print("="*70)
    
    print(f"\nM-CSA Ground Truth: {sorted(true_set)}")
    print(f"Predicted (Top {len(predicted)}): {sorted(predicted)}")
    
    print("\n--- Detailed Analysis ---")
    
    # True Positives
    tp_positions = predicted_set & true_set
    if tp_positions:
        print(f"\n✓ TRUE POSITIVES ({len(tp_positions)}):")
        for pos in sorted(tp_positions):
            pos_data = _get_pos_data(pos)
            if pos_data:
                print(f"   Residue {pos} ({pos_data['consensus']}): "
                      f"conservation={pos_data['conservation']:.4f}")
    
    # False Positives
    fp_positions = predicted_set - true_set
    if fp_positions:
        print(f"\n✗ FALSE POSITIVES ({len(fp_positions)}):")
        for pos in sorted(fp_positions):
            pos_data = _get_pos_data(pos)
            if pos_data:
                print(f"   Residue {pos} ({pos_data['consensus']}): "
                      f"conservation={pos_data['conservation']:.4f}")
    
    # False Negatives
    fn_positions = true_set - predicted_set
    if fn_positions:
        print(f"\n✗ FALSE NEGATIVES (Missed catalytic residues: {len(fn_positions)}):")
        for pos in sorted(fn_positions):
            pos_data = _get_pos_data(pos)
            if pos_data:
                rank = get_conservation_rank(conservation_data, pos, alignment_mapping)
                print(f"   Residue {pos} ({pos_data['consensus']}): "
                      f"conservation={pos_data['conservation']:.4f}, "
                      f"rank={rank}")
            else:
                print(f"   Residue {pos}: not found in conservation data")
    
    print("\n--- Performance Metrics ---")
    print(f"Precision: {metrics['precision']:.3f} ({metrics['tp']}/{metrics['n_predicted']})")
    print(f"Recall:    {metrics['recall']:.3f} ({metrics['tp']}/{metrics['n_true']})")
    print(f"F1 Score:  {metrics['f1']:.3f}")
    
    print("\n" + "="*70)
    print("This is your BASELINE. More complex models must beat this.")
    print("="*70)

def get_conservation_rank(conservation_data: Dict, resid: int, 
                          alignment_mapping: Dict = None) -> int:
    """Get the rank of a residue (auth_resid) by conservation score."""
    # Build resid -> conservation score mapping
    if alignment_mapping:
        col_to_resid = {int(k): v for k, v in alignment_mapping.items()}
    else:
        col_to_resid = {}
    
    # Pair each position with its resid
    scored = []
    for p in conservation_data['positions']:
        r = p.get('resid') or p.get('auth_resid') or col_to_resid.get(p['position'])
        if r is not None:
            scored.append((r, p['conservation']))
    
    # Sort by conservation descending
    scored.sort(key=lambda x: x[1], reverse=True)
    
    for rank, (r, _) in enumerate(scored, 1):
        if r == resid:
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
    parser.add_argument('--top-n', default='auto',
                        help='Top N predictions. "auto" matches ground truth count (default: auto)')
    parser.add_argument('--conservation-threshold', type=float, default=None,
                        help='Conservation score threshold (overrides --top-n). '
                             'E.g., 0.6 selects all residues with conservation >= 0.6')
    parser.add_argument('--exclude-gaps', action='store_true', default=True)
    parser.add_argument('--exclude-structural', action='store_true', default=False,
                        help='Hard-filter G, P, A residues (not recommended, use --catalytic-propensity instead)')
    parser.add_argument('--catalytic-propensity', action='store_true', default=False,
                        help='Weight conservation by M-CSA catalytic propensity (recommended). '
                             'Boosts H, C, D, E; dampens G, P, A, hydrophobics.')
    parser.add_argument('--p2rank-json', default=None,
                        help='P2Rank scores JSON from parse_p2rank.py (binding site prediction)')
    parser.add_argument('--pdb-file', default=None,
                        help='Path to PDB file for spatial clustering (optional)')
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
    
    # Resolve top-n: auto = match ground truth count
    if args.conservation_threshold is not None:
        top_n = None  # threshold mode
        print(f"  Using conservation threshold: {args.conservation_threshold}")
    elif args.top_n == 'auto':
        top_n = len(true_positions)
        print(f"  Auto top-N: predicting {top_n} residues (matching ground truth count)")
    else:
        top_n = int(args.top_n)
    
    # Load P2Rank scores if provided
    p2rank_scores = None
    if args.p2rank_json:
        print(f"Loading P2Rank scores: {args.p2rank_json}")
        with open(args.p2rank_json) as f:
            p2rank_data = json.load(f)
            p2rank_scores = p2rank_data.get('residues', {})
        print(f"  Loaded scores for {len(p2rank_scores)} residues")
    
    # Get predictions
    predicted_positions = get_top_conserved_positions(
        conservation_data,
        alignment_mapping,
        top_n=top_n if top_n else len(true_positions),
        conservation_threshold=args.conservation_threshold,
        exclude_gaps=args.exclude_gaps,
        min_identity=args.min_identity,
        exclude_structural=args.exclude_structural,
        use_catalytic_propensity=args.catalytic_propensity,
        p2rank_scores=p2rank_scores,
        pdb_file=args.pdb_file
    )
    
    # Calculate metrics
    metrics = calculate_metrics(predicted_positions, true_positions)
    
    # Print results
    print_results(predicted_positions, true_positions, metrics, conservation_data, alignment_mapping)
    
    # Save to file if requested
    if args.output:
        results = {
            'pdb_id': args.pdb_id,
            'top_n': top_n,
            'conservation_threshold': args.conservation_threshold,
            'mcsa_ground_truth': sorted(true_positions),
            'predicted': sorted(predicted_positions),
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