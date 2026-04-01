# Download API - Catalytic residues
# - Each array element = one catalytic residue. Residues belonging to the same enzyme share the same mcsa_id.
wget https://www.ebi.ac.uk/thornton-srv/m-csa/api/residues/?format=json -O m-csa_residues.json

{
"mcsa_id": 2,
"roles_summary": "activator, electrostatic stabiliser, hydrogen bond acceptor, hydrogen bond donor, proton acceptor, proton donor",
"function_location_abv": "",
"main_annotation": "In acylation the residue acts as a general base towards a structurally conserved water molecule,leading to the deprotonation of Ser70 (proton relay). In deacylation, Glu 166 abstracts a proton from a water molecule, activating a nucleophile for attack at the substrate carbon linked to the gamma oxygen of S 70. The residue is hydrogen bonded to Lys 73.",
"ptm": "",
"roles": [
    {
    "group_function": "",
    "function_type": "interaction",
    "function": "hydrogen bond acceptor",
    "emo": "EMO_00113"
    },
    {
    "group_function": "activator",
    "function_type": "spectator",
    "function": "activator",
    "emo": "EMO_00038"
    },
    {
    "group_function": "",
    "function_type": "interaction",
    "function": "hydrogen bond acceptor",
    "emo": "EMO_00113"
    },
    {
    "group_function": "electrostatic interaction",
    "function_type": "spectator",
    "function": "electrostatic stabiliser",
    "emo": "EMO_00033"
    },
    {
    "group_function": "",
    "function_type": "interaction",
    "function": "hydrogen bond acceptor",
    "emo": "EMO_00113"
    },
    {
    "group_function": "",
    "function_type": "interaction",
    "function": "hydrogen bond acceptor",
    "emo": "EMO_00113"
    },
    {
    "group_function": "",
    "function_type": "interaction",
    "function": "hydrogen bond donor",
    "emo": "EMO_00114"
    },
    {
    "group_function": "",
    "function_type": "interaction",
    "function": "hydrogen bond acceptor",
    "emo": "EMO_00113"
    },
    {
    "group_function": "",
    "function_type": "interaction",
    "function": "hydrogen bond acceptor",
    "emo": "EMO_00113"
    },
    {
    "group_function": "activator",
    "function_type": "spectator",
    "function": "activator",
    "emo": "EMO_00038"
    },
    {
    "group_function": "",
    "function_type": "interaction",
    "function": "hydrogen bond acceptor",
    "emo": "EMO_00113"
    },
    {
    "group_function": "electrostatic interaction",
    "function_type": "spectator",
    "function": "electrostatic stabiliser",
    "emo": "EMO_00033"
    },
    {
    "group_function": "",
    "function_type": "interaction",
    "function": "hydrogen bond acceptor",
    "emo": "EMO_00113"
    },
    {
    "group_function": "",
    "function_type": "interaction",
    "function": "hydrogen bond acceptor",
    "emo": "EMO_00113"
    },
    {
    "group_function": "",
    "function_type": "interaction",
    "function": "hydrogen bond donor",
    "emo": "EMO_00114"
    },
    {
    "group_function": "",
    "function_type": "interaction",
    "function": "hydrogen bond acceptor",
    "emo": "EMO_00113"
    },
    {
    "group_function": "proton shuttle (general acid/base)",
    "function_type": "reactant",
    "function": "proton acceptor",
    "emo": "EMO_00066"
    },
    {
    "group_function": "proton shuttle (general acid/base)",
    "function_type": "reactant",
    "function": "proton acceptor",
    "emo": "EMO_00066"
    },
    {
    "group_function": "proton shuttle (general acid/base)",
    "function_type": "reactant",
    "function": "proton acceptor",
    "emo": "EMO_00066"
    },
    {
    "group_function": "proton shuttle (general acid/base)",
    "function_type": "reactant",
    "function": "proton donor",
    "emo": "EMO_00068"
    },
    {
    "group_function": "proton shuttle (general acid/base)",
    "function_type": "reactant",
    "function": "proton donor",
    "emo": "EMO_00068"
    },
    {
    "group_function": "proton shuttle (general acid/base)",
    "function_type": "reactant",
    "function": "proton donor",
    "emo": "EMO_00068"
    }
],
"residue_chains": [ # the residue in PDB structure space
    {
    "chain_name": "A", # chain in the mmCIF file
    "pdb_id": "1btl", # which PDB structure
    "assembly_chain_name": "A", # chain in biological assembly
    "assembly": 1, # biological assembly number
    "code": "Glu", # amino acid
    "resid": 141, # residue number in mmCIF file
    "auth_resid": 166, # residue number in PDB file - what FoldDisco's -q flag expects??
    "is_reference": true, # is this the reference PDB for this enzyme
    "domain_name": "A00", # CATH domain name
    "domain_cath_id": "3.40.710.10" # CATH domain ID
    }
],
"residue_sequences": [ # the same residue in UniProt sequence space
    {
    "uniprot_id": "P62593", # UniProt accession
    "code": "Glu",
    "is_reference": true,
    "resid": 164 # position in the UniProt sequence
    }
]
},