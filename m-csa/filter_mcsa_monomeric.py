"""
Filter M-CSA database entries for monomeric, single-chain protein structures
by querying the RCSB PDB Search API.

Usage:
    python filter_mcsa_monomeric.py --input mcsa_parsed.tsv --output mcsa_monomeric.tsv

Optionally validate a random sample of results:
    python filter_mcsa_monomeric.py --input mcsa_parsed.tsv --output mcsa_monomeric.tsv --validate
"""

import argparse
import csv
import json
import random
import requests


SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"


def get_monomeric_protein_pdb_ids() -> set:
    """
    Query RCSB Search API for all PDB entries where:
      1. Exactly 1 protein entity (one unique protein sequence)
      2. Exactly 1 polymer chain deposited in the coordinates
      3. No DNA or RNA entities

    We use entry-level attributes only. `rcsb_struct_symmetry.oligomeric_state`
    would be stricter (biological-assembly monomer) but is a nested-array field
    that returns empty results via simple operators — for M-CSA reference
    structures the deposited coords usually match the biological assembly, so
    entry-level filtering is a good-enough proxy.

    Returns a set of uppercase PDB IDs.
    """
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
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
                        "attribute": "rcsb_entry_info.deposited_polymer_entity_instance_count",
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

    # RCSB returns 204 No Content when the query is valid but matches zero entries.
    if resp.status_code == 204:
        print("WARNING: query returned zero hits — check attribute names/values")
        return set()
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

    in_sep = "\t" if input_path.endswith(".tsv") else ","
    with open(input_path, newline="") as f:
        reader = csv.DictReader(f, delimiter=in_sep)
        fieldnames = reader.fieldnames
        rows = list(reader)
    print(f"Loaded {len(rows)} rows from {input_path}")

    mcsa_pdb_ids = {r["pdb_id"].upper() for r in rows if r.get("pdb_id")}
    print(f"Unique PDB IDs in M-CSA: {len(mcsa_pdb_ids)}")

    # Get monomeric protein PDB IDs from RCSB
    mono_ids = get_monomeric_protein_pdb_ids()

    # Intersect
    included = mcsa_pdb_ids & mono_ids
    excluded = mcsa_pdb_ids - mono_ids

    print(f"\nResults:")
    print(f"  Monomeric protein entries in M-CSA: {len(included)}")
    print(f"  Excluded: {len(excluded)}")

    # Filter rows and save
    filtered = [r for r in rows if r.get("pdb_id", "").upper() in mono_ids]
    print(f"  Rows in filtered output: {len(filtered)}")

    out_sep = "\t" if output_path.endswith(".tsv") else ","
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter=out_sep)
        writer.writeheader()
        writer.writerows(filtered)
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