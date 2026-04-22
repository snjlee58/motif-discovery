"""
Filter M-CSA database entries for biological monomers (exactly 1 polymer chain
in the biological assembly, a single protein entity, no DNA/RNA) by querying
the RCSB PDB Data API per-entry.

Why per-entry instead of one Search API call? The Search API's entry-level
attributes can't distinguish deposited-ASU chain count from biological-assembly
chain count. Homodimers crystallized with 1 chain/ASU would pass a
deposited-count filter, and true monomers deposited with NCS copies would fail
it. The Data API's /assembly/{pdb}/1 endpoint reports the correct biological
chain count, so we verify each M-CSA entry directly.

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
from concurrent.futures import ThreadPoolExecutor, as_completed


ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry"
ASSEMBLY_URL = "https://data.rcsb.org/rest/v1/core/assembly"


def is_biological_monomer(pdb_id: str) -> bool:
    """
    Return True iff `pdb_id` meets all of:
      - biological assembly (assembly 1) has exactly 1 polymer chain instance
      - exactly 1 protein entity
      - no DNA or RNA entities

    Makes two Data API calls (entry + assembly/1). Returns False on any failure
    so transient network errors don't corrupt the output — rerun to retry.
    """
    session = _session()
    try:
        entry_resp = session.get(f"{ENTRY_URL}/{pdb_id}", timeout=30)
        if entry_resp.status_code != 200:
            return False
        info = entry_resp.json().get("rcsb_entry_info", {})
        if info.get("polymer_entity_count_protein") != 1:
            return False
        if info.get("polymer_entity_count_DNA", 0) != 0:
            return False
        if info.get("polymer_entity_count_RNA", 0) != 0:
            return False

        asm_resp = session.get(f"{ASSEMBLY_URL}/{pdb_id}/1", timeout=30)
        if asm_resp.status_code != 200:
            return False
        asm_info = asm_resp.json().get("rcsb_assembly_info", {})
        return asm_info.get("polymer_entity_instance_count") == 1
    except requests.RequestException:
        return False


# One session per thread — lets urllib3 reuse HTTP connections.
import threading
_thread_local = threading.local()
def _session() -> requests.Session:
    s = getattr(_thread_local, "session", None)
    if s is None:
        s = requests.Session()
        _thread_local.session = s
    return s


def classify_pdb_ids(pdb_ids: set, max_workers: int = 20) -> tuple:
    """
    Check each PDB ID against the Data API in parallel.

    Returns:
        (included, excluded) — both sets of uppercase PDB IDs.
    """
    pdb_list = sorted(pdb_ids)
    included = set()
    total = len(pdb_list)
    done = 0
    progress_every = max(1, total // 20)

    print(f"Verifying biological monomer status for {total} entries "
          f"(Data API, {max_workers} threads)...")

    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = {ex.submit(is_biological_monomer, pid): pid for pid in pdb_list}
        for fut in as_completed(futures):
            pid = futures[fut]
            try:
                if fut.result():
                    included.add(pid)
            except Exception:
                pass  # leave in excluded
            done += 1
            if done % progress_every == 0 or done == total:
                print(f"  {done}/{total} checked ({len(included)} monomers so far)")

    excluded = pdb_ids - included
    return included, excluded


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

    # Per-entry Data API check (biological assembly chain count)
    included, excluded = classify_pdb_ids(mcsa_pdb_ids)

    print(f"\nResults:")
    print(f"  Biological monomers in M-CSA: {len(included)}")
    print(f"  Excluded: {len(excluded)}")

    # Filter rows and save
    filtered = [r for r in rows if r.get("pdb_id", "").upper() in included]
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