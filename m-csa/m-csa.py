import json

with open("/mnt/scratch/sunny/m-csa/mcsa_residues.json") as f:
    data = json.load(f)

pdb_ids = set()
uniprot_ids = set()

for residue in data:
    for chain in residue.get('residue_chains', []):
        if chain.get('is_reference'):
            pdb_ids.add(chain['pdb_id'].lower())
    for seq in residue.get('residue_sequences', []):
        if seq.get('is_reference'):
            uniprot_ids.add(seq['uniprot_id'])

print(f"Total residue entries:     {len(data)}")
print(f"Unique reference PDB IDs:  {len(pdb_ids)}")
print(f"Unique reference UniProts: {len(uniprot_ids)}")