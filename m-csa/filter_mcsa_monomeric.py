"""
Filter M-CSA database entries for monomeric, single-chain protein structures
by querying the RCSB PDB Search API.

Usage:
    python filter_mcsa_monomeric.py --input mcsa_parsed.tsv --output mcsa_monomeric.tsv

Optionally validate a random sample of results:
    python filter_mcsa_monomeric.py --input mcsa_parsed.tsv --output mcsa_monomeric.tsv --validate
"""

import argparse
import json
import random
import requests
import pandas as pd


SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"


def get_monomeric_protein_pdb_ids() -> set:
    """
    Query RCSB Search API for all PDB entries where:
      1. RCSB's curated biological-assembly oligomeric state is "Monomeric"
      2. The entry has exactly 1 protein entity
      3. The entry has no DNA or RNA entities

    Returns a set of uppercase PDB IDs.
    """
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    # Use the curated oligomeric state rather than raw chain-instance count.
                    # rcsb_assembly_info.polymer_entity_instance_count is per-assembly and
                    # matches the deposited ASU (not just the biological assembly), so it
                    # incorrectly passes homodimers that crystallize with one chain per ASU.
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_struct_symmetry.oligomeric_state",
                        "operator": "equals",
                        "value": "Monomeric"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.polymer_entity_count_protein",
                        "operator": "equals",
                        "value": 1
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.polymer_entity_count_DNA",
                        "operator": "equals",
                        "value": 0
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.polymer_entity_count_RNA",
                        "operator": "equals",
                        "value": 0
                    }
                }
            ]
        },
        "return_type": "entry",
        "request_options": {
            "return_all_hits": True
        }
    }

    print("Querying RCSB Search API for monomeric protein entries...")
    resp = requests.post(SEARCH_URL, json=query, timeout=120)
    print(f"HTTP {resp.status_code}")

    if resp.status_code != 200:
        print(f"Response: {resp.text[:500]}")
        raise RuntimeError(f"RCSB Search API failed with HTTP {resp.status_code}")

    data = resp.json()
    total = data.get("total_count", 0)
    print(f"Total monomeric protein entries: {total}")

    # return_type is "entry", so identifiers are plain PDB IDs like "105M"
    pdb_ids = set()
    for hit in data.get("result_set", []):
        pdb_ids.add(hit["identifier"].upper())

    return pdb_ids


def validate_results(included: set, excluded: set, n: int = 10):
    """
    Spot-check a random sample of included and excluded PDB IDs
    by querying the RCSB Data API for their actual chain/entity counts.
    """
    print(f"\n{'='*60}")
    print("VALIDATION: Spot-checking results against RCSB Data API")
    print(f"{'='*60}")

    included_sample = random.sample(sorted(included), min(n, len(included)))
    excluded_sample = random.sample(sorted(excluded), min(n, len(excluded)))

    def check_entry(pdb_id):
        # Entry-level info
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        resp = requests.get(url, timeout=30)
        if resp.status_code != 200:
            print(f"  {pdb_id}: could not fetch entry data (HTTP {resp.status_code})")
            return

        info = resp.json().get("rcsb_entry_info", {})

        # Assembly-level info (assembly 1)
        asm_url = f"https://data.rcsb.org/rest/v1/core/assembly/{pdb_id}/1"
        asm_resp = requests.get(asm_url, timeout=30)
        asm_chains = "?"
        if asm_resp.status_code == 200:
            asm_info = asm_resp.json().get("rcsb_assembly_info", {})
            asm_chains = asm_info.get("polymer_entity_instance_count", "?")

        protein = info.get("polymer_entity_count_protein", "?")
        dna = info.get("polymer_entity_count_DNA", "?")
        rna = info.get("polymer_entity_count_RNA", "?")

        print(f"  {pdb_id}: chains={asm_chains}  protein={protein}  DNA={dna}  RNA={rna}")

    print(f"\nINCLUDED (expect: chains=1, protein=1, DNA=0, RNA=0):")
    for pid in included_sample:
        check_entry(pid)

    print(f"\nEXCLUDED (expect: at least one condition violated):")
    for pid in excluded_sample:
        check_entry(pid)


def filter_mcsa(input_path: str, output_path: str, validate: bool = False):
    """Read M-CSA parsed data, filter for monomeric entries, and save."""

    sep = "\t" if input_path.endswith(".tsv") else ","
    df = pd.read_csv(input_path, sep=sep)
    print(f"Loaded {len(df)} rows from {input_path}")

    mcsa_pdb_ids = set(df["pdb_id"].str.upper().unique())
    print(f"Unique PDB IDs in M-CSA: {len(mcsa_pdb_ids)}")

    # Get monomeric protein PDB IDs from RCSB
    mono_ids = get_monomeric_protein_pdb_ids()

    # Intersect
    included = mcsa_pdb_ids & mono_ids
    excluded = mcsa_pdb_ids - mono_ids

    print(f"\nResults:")
    print(f"  Monomeric protein entries in M-CSA: {len(included)}")
    print(f"  Excluded: {len(excluded)}")

    # Filter and save
    df_filtered = df[df["pdb_id"].str.upper().isin(mono_ids)].copy()
    print(f"  Rows in filtered dataframe: {len(df_filtered)}")

    out_sep = "\t" if output_path.endswith(".tsv") else ","
    df_filtered.to_csv(output_path, index=False, sep=out_sep)
    print(f"\nSaved to {output_path}")

    # Save PDB ID lists
    base = output_path.rsplit(".", 1)[0]

    with open(f"{base}_pdb_ids.txt", "w") as f:
        for pid in sorted(included):
            f.write(pid + "\n")
    print(f"Included PDB IDs -> {base}_pdb_ids.txt")

    with open(f"{base}_excluded_pdb_ids.txt", "w") as f:
        for pid in sorted(excluded):
            f.write(pid + "\n")
    print(f"Excluded PDB IDs -> {base}_excluded_pdb_ids.txt")

    # Show examples
    if included:
        print(f"\nExample included: {sorted(included)[:10]}")
    if excluded:
        print(f"Example excluded: {sorted(excluded)[:10]}")

    # Optional validation
    if validate:
        validate_results(included, excluded)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter M-CSA entries for monomeric single-chain proteins"
    )
    parser.add_argument("--input", required=True, help="Path to parsed M-CSA TSV/CSV")
    parser.add_argument("--output", required=True, help="Path for filtered output")
    parser.add_argument("--validate", action="store_true",
                        help="Spot-check 10 included + 10 excluded IDs against Data API")
    args = parser.parse_args()

    filter_mcsa(args.input, args.output, args.validate)