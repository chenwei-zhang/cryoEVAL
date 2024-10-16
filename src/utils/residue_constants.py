######################################################
### Credit from ModelAngelo evaluation functions  ####
######################################################


"""
Vast majority of this file comes from 
http://github.com/deepmind/alphafold/blob/main/alphafold/common/residue_constants.py

Nucleotides come from
https://github.com/uw-ipd/RoseTTAFold2NA/blob/main/network/chemical.py
"""
import functools
import os
from collections import namedtuple
from typing import List, Mapping, Tuple

import numpy as np
import torch
import einops

restype_1to3 = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
    "x": "DA",
    "y": "DC",
    "z": "DG",
    "t": "DT",
    "a": "A",
    "c": "C",
    "g": "G",
    "u": "U",
    "1": "UNK",
    "2": "DN",
    "3": "N",
}

restype_3to1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "DA": "x",
    "DC": "y",
    "DG": "z",
    "DT": "t",
    "A": "a",
    "C": "c",
    "G": "g",
    "U": "u",
    "UNK": "1",
    "DN": "2",
    "N": "3",
}

restype_3_to_index = {
    "ALA": 0,
    "ARG": 1,
    "ASN": 2,
    "ASP": 3,
    "CYS": 4,
    "GLN": 5,
    "GLU": 6,
    "GLY": 7,
    "HIS": 8,
    "ILE": 9,
    "LEU": 10,
    "LYS": 11,
    "MET": 12,
    "PHE": 13,
    "PRO": 14,
    "SER": 15,
    "THR": 16,
    "TRP": 17,
    "TYR": 18,
    "VAL": 19,
    "DA": 20,
    "DC": 21,
    "DG": 22,
    "DT": 23,
    "A": 24,
    "C": 25,
    "G": 26,
    "U": 27,
    "UNK":28,
    "DN": 29,
    "N": 30,
}

restype_1_to_index = {
    "A": 0,
    "R": 1,
    "N": 2,
    "D": 3,
    "C": 4,
    "Q": 5,
    "E": 6,
    "G": 7,
    "H": 8,
    "I": 9,
    "L": 10,
    "K": 11,
    "M": 12,
    "F": 13,
    "P": 14,
    "S": 15,
    "T": 16,
    "W": 17,
    "Y": 18,
    "V": 19,
    "x": 20,
    "y": 21,
    "z": 22,
    "t": 23,
    "a": 24,
    "c": 25,
    "g": 26,
    "u": 27,
    "1": 28,
    "2": 29,
    "3": 30,
}

index_to_restype_1 = [
    "A",
    "R",
    "N",
    "D",
    "C",
    "Q",
    "E",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V",
    "x",
    "y",
    "z",
    "t",
    "a",
    "c",
    "g",
    "u",
    "1",
    "2",
    "3",
]

index_to_restype_3 = [
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "DA",
    "DC",
    "DG",
    "DT",
    "A",
    "C",
    "G",
    "U",
    "UNK",
    "DN",
    "N",
]

num_prot = 20
prot_restype3 = set(index_to_restype_3[:num_prot] + ["UNK"])
prot_restype1 = set(index_to_restype_1[:num_prot] + ["1"])

index_to_nuc = index_to_restype_3[num_prot:]
index_to_nuc.pop(index_to_nuc.index("UNK"))


def restype3_is_na(restype3: str) -> bool:
    return not (restype3 in prot_restype3)


def restype1_is_na(restype1: str) -> bool:
    return not (restype1 in prot_restype1)


def restype3_is_prot(restype3: str) -> bool:
    return restype3 in prot_restype3


def restype1_is_prot(restype1: str) -> bool:
    return restype1 in prot_restype1


index_to_hmm_restype_1 = sorted(index_to_restype_1[:num_prot])
hmm_restype_1_to_index = {
    hmm_restype: i for i, hmm_restype in enumerate(index_to_hmm_restype_1)
}
restype_1_order_to_hmm = [restype_1_to_index[aa] for aa in index_to_hmm_restype_1]

restype_order = {restype: i for i, restype in enumerate(index_to_restype_1)}
restype_num = len(index_to_restype_1)

unk_restype_index = restype_num  # Catch-all index for unknown restypes.

restypes_with_x = index_to_restype_1 + ["X"]
restype_order_with_x = {restype: i for i, restype in enumerate(restypes_with_x)}

# This mapping is used when we need to store atom data in a format that requires
# fixed atom data size for every residue (e.g. a numpy array).
atom_types = [
    "N",
    "CA",
    "C",
    "O",
    "CB",
    "CG",
    "CD",
    "NE",
    "CZ",
    "NH1",
    "NH2",
    "OD1",
    "ND2",
    "OD2",
    "SG",
    "OE1",
    "NE2",
    "OE2",
    "ND1",
    "CD2",
    "CE1",
    "CG1",
    "CG2",
    "CD1",
    "CE",
    "NZ",
    "SD",
    "CE2",
    "OG",
    "OG1",
    "NE1",
    "CE3",
    "CZ2",
    "CZ3",
    "CH2",
    "OH",
    "OP1",
    "P",
    "OP2",
    "O5'",
    "C5'",
    "C4'",
    "O4'",
    "C3'",
    "O3'",
    "C2'",
    "C1'",
    "N9",
    "C4",
    "N3",
    "C2",
    "N1",
    "C6",
    "C5",
    "N7",
    "C8",
    "N6",
    "O2",
    "N4",
    "N2",
    "O6",
    "O4",
    "C7",
    "O2'",
    "OXT",
]

atom_order = {atom_type: i for i, atom_type in enumerate(atom_types)}
atom_type_num = len(atom_types)  # := 65.
num_atoms = atom_type_num

# atomf -> atom full
atomf_backbone_mask = np.zeros((1, atom_type_num), dtype=np.float32)
atomf_backbone_mask[
    :,
    [
        atom_order["N"],
        atom_order["CA"],
        atom_order["C"],
        atom_order["O"],
        atom_order["OP1"],
        atom_order["P"],
        atom_order["OP2"],
        atom_order["O5'"],
    ],
] = 1

cif_secondary_structure_to_index = {
    "NULL": 0,
    "HELIX_1": 1,
    "HELIX_5": 2,
    "SHEET": 3,
}

ca_to_n_distance_ang = 1.4556349
ca_to_c_distance_ang = 1.5235157
peptide_bond_length_ang = 1.3310018

c_to_ca_to_n_angle_rad = 1.9384360


def parse_sequence_string(sequence: str) -> List[int]:
    return [restype_1_to_index[s] for s in sequence]


def parse_index_list(index_list: List[int]) -> List[str]:
    return "".join(index_to_restype_1[i] for i in index_list)


def translate_restype_3_to_1(residue_list: List[str]) -> List[str]:
    return [restype_3to1[r] for r in residue_list]


restype3_to_atoms = {
    "ALA": ["N", "CA", "C", "O", "CB"],
    "ARG": ["N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
    "ASN": ["N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"],
    "ASP": ["N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"],
    "CYS": ["N", "CA", "C", "O", "CB", "SG"],
    "GLN": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"],
    "GLU": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"],
    "GLY": ["N", "CA", "C", "O"],
    "HIS": ["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"],
    "ILE": ["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"],
    "LEU": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"],
    "LYS": ["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"],
    "MET": ["N", "CA", "C", "O", "CB", "CG", "SD", "CE"],
    "PHE": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "PRO": ["N", "CA", "C", "O", "CB", "CG", "CD"],
    "SER": ["N", "CA", "C", "O", "CB", "OG"],
    "THR": ["N", "CA", "C", "O", "CB", "OG1", "CG2"],
    "TRP": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "CD1",
        "CD2",
        "NE1",
        "CE2",
        "CE3",
        "CZ2",
        "CZ3",
        "CH2",
    ],
    "TYR": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
    "VAL": ["N", "CA", "C", "O", "CB", "CG1", "CG2"],
    "DA": [
        "OP1",
        "P",
        "OP2",
        "O5'",
        "C5'",
        "C4'",
        "O4'",
        "C3'",
        "O3'",
        "C2'",
        "C1'",
        "N9",
        "C4",
        "N3",
        "C2",
        "N1",
        "C6",
        "C5",
        "N7",
        "C8",
        "N6",
    ],
    "DC": [
        "OP1",
        "P",
        "OP2",
        "O5'",
        "C5'",
        "C4'",
        "O4'",
        "C3'",
        "O3'",
        "C2'",
        "C1'",
        "N1",
        "C2",
        "O2",
        "N3",
        "C4",
        "N4",
        "C5",
        "C6",
    ],
    "DG": [
        "OP1",
        "P",
        "OP2",
        "O5'",
        "C5'",
        "C4'",
        "O4'",
        "C3'",
        "O3'",
        "C2'",
        "C1'",
        "N9",
        "C4",
        "N3",
        "C2",
        "N1",
        "C6",
        "C5",
        "N7",
        "C8",
        "N2",
        "O6",
    ],
    "DT": [
        "OP1",
        "P",
        "OP2",
        "O5'",
        "C5'",
        "C4'",
        "O4'",
        "C3'",
        "O3'",
        "C2'",
        "C1'",
        "N1",
        "C2",
        "O2",
        "N3",
        "C4",
        "O4",
        "C5",
        "C7",
        "C6",
    ],
    "A": [
        "OP1",
        "P",
        "OP2",
        "O5'",
        "C5'",
        "C4'",
        "O4'",
        "C3'",
        "O3'",
        "C1'",
        "C2'",
        "O2'",
        "N1",
        "C2",
        "N3",
        "C4",
        "C5",
        "C6",
        "N6",
        "N7",
        "C8",
        "N9",
    ],
    "C": [
        "OP1",
        "P",
        "OP2",
        "O5'",
        "C5'",
        "C4'",
        "O4'",
        "C3'",
        "O3'",
        "C1'",
        "C2'",
        "O2'",
        "N1",
        "C2",
        "O2",
        "N3",
        "C4",
        "N4",
        "C5",
        "C6",
    ],
    "G": [
        "OP1",
        "P",
        "OP2",
        "O5'",
        "C5'",
        "C4'",
        "O4'",
        "C3'",
        "O3'",
        "C1'",
        "C2'",
        "O2'",
        "N1",
        "C2",
        "N2",
        "N3",
        "C4",
        "C5",
        "C6",
        "O6",
        "N7",
        "C8",
        "N9",
    ],
    "U": [
        "OP1",
        "P",
        "OP2",
        "O5'",
        "C5'",
        "C4'",
        "O4'",
        "C3'",
        "O3'",
        "C1'",
        "C2'",
        "O2'",
        "N1",
        "C2",
        "O2",
        "N3",
        "C4",
        "O4",
        "C5",
        "C6",
    ],
}
restype3_to_atoms["UNK"] = restype3_to_atoms["GLY"]
restype3_to_atoms["DN"] = restype3_to_atoms["DC"]
restype3_to_atoms["N"] = restype3_to_atoms["C"]

restype3_to_atoms_index = dict(
    [
        (res, dict([(a, i) for (i, a) in enumerate(atoms)]))
        for (res, atoms) in restype3_to_atoms.items()
    ]
)
for residue in restype3_to_atoms_index:
    if restype3_is_prot(residue):
        restype3_to_atoms_index[residue]["OXT"] = restype3_to_atoms_index[residue]["O"]

full_num_residues = len(restype3_to_atoms_index)
canonical_num_residues = full_num_residues - 3
backbone_atoms_prot = {"CA", "C", "N"}
backbone_atoms_nuc = {"OP1", "P", "OP2"}
backbone_atoms = backbone_atoms_prot.union(backbone_atoms_nuc)

secondary_structure_to_simplified_index = {
    # CIF
    "NULL": 0,
    "HELIX_1": 1,
    "HELIX_5": 1,
    "SHEET": 2,
    # DSSP
    "OTHER": 0,
    "BEND": 0,
    "TURN_TY1_P": 0,
    "HELX_RH_AL_P": 1,
    "HELX_RH_3T_P": 1,
    "HELX_LH_PP_P": 1,
    "HELX_RH_PI_P": 1,
    "STRN": 2,
}

# Distance from one CA to next CA [trans configuration: omega = 180].
ca_ca = 3.80209737096

# Format: The list for each AA type contains chi1, chi2, chi3, chi4 in
# this order (or a relevant subset from chi1 onwards). ALA and GLY don't have
# chi angles so their chi angle lists are empty.
chi_angles_atoms = {
    "ALA": [],
    # Chi5 in arginine is always 0 +- 5 degrees, so ignore it.
    "ARG": [
        ["N", "CA", "CB", "CG"],
        ["CA", "CB", "CG", "CD"],
        ["CB", "CG", "CD", "NE"],
        ["CG", "CD", "NE", "CZ"],
    ],
    "ASN": [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "OD1"]],
    "ASP": [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "OD1"]],
    "CYS": [["N", "CA", "CB", "SG"]],
    "GLN": [
        ["N", "CA", "CB", "CG"],
        ["CA", "CB", "CG", "CD"],
        ["CB", "CG", "CD", "OE1"],
    ],
    "GLU": [
        ["N", "CA", "CB", "CG"],
        ["CA", "CB", "CG", "CD"],
        ["CB", "CG", "CD", "OE1"],
    ],
    "GLY": [],
    "HIS": [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "ND1"]],
    "ILE": [["N", "CA", "CB", "CG1"], ["CA", "CB", "CG1", "CD1"]],
    "LEU": [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]],
    "LYS": [
        ["N", "CA", "CB", "CG"],
        ["CA", "CB", "CG", "CD"],
        ["CB", "CG", "CD", "CE"],
        ["CG", "CD", "CE", "NZ"],
    ],
    "MET": [
        ["N", "CA", "CB", "CG"],
        ["CA", "CB", "CG", "SD"],
        ["CB", "CG", "SD", "CE"],
    ],
    "PHE": [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]],
    "PRO": [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"]],
    "SER": [["N", "CA", "CB", "OG"]],
    "THR": [["N", "CA", "CB", "OG1"]],
    "TRP": [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]],
    "TYR": [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]],
    "VAL": [["N", "CA", "CB", "CG1"]],
    "DA": [["C3'", "C2'", "N9", "C4"]],
    "DC": [["C3'", "C2'", "N1", "C2"]],
    "DG": [["C3'", "C2'", "N9", "C4"]],
    "DT": [["C3'", "C2'", "N1", "C2"]],
    "A": [["C3'", "C2'", "N9", "C4"]],
    "C": [["C3'", "C2'", "N1", "C2"]],
    "G": [["C3'", "C2'", "N9", "C4"]],
    "U": [["C3'", "C2'", "N1", "C2"]],
    "UNK": [],
    "DN": [["C3'", "C2'", "N1", "C2"]],
    "N": [["C3'", "C2'", "N1", "C2"]],
}

# If chi angles given in fixed-length array, this matrix determines how to mask
# them for each AA type. The order is as per restype_order (see below).
chi_angles_mask = [
    [0.0, 0.0, 0.0, 0.0],  # ALA
    [1.0, 1.0, 1.0, 1.0],  # ARG
    [1.0, 1.0, 0.0, 0.0],  # ASN
    [1.0, 1.0, 0.0, 0.0],  # ASP
    [1.0, 0.0, 0.0, 0.0],  # CYS
    [1.0, 1.0, 1.0, 0.0],  # GLN
    [1.0, 1.0, 1.0, 0.0],  # GLU
    [0.0, 0.0, 0.0, 0.0],  # GLY
    [1.0, 1.0, 0.0, 0.0],  # HIS
    [1.0, 1.0, 0.0, 0.0],  # ILE
    [1.0, 1.0, 0.0, 0.0],  # LEU
    [1.0, 1.0, 1.0, 1.0],  # LYS
    [1.0, 1.0, 1.0, 0.0],  # MET
    [1.0, 1.0, 0.0, 0.0],  # PHE
    [1.0, 1.0, 0.0, 0.0],  # PRO
    [1.0, 0.0, 0.0, 0.0],  # SER
    [1.0, 0.0, 0.0, 0.0],  # THR
    [1.0, 1.0, 0.0, 0.0],  # TRP
    [1.0, 1.0, 0.0, 0.0],  # TYR
    [1.0, 0.0, 0.0, 0.0],  # VAL
    [1.0, 0.0, 0.0, 0.0],  # DA
    [1.0, 0.0, 0.0, 0.0],  # DC
    [1.0, 0.0, 0.0, 0.0],  # DG
    [1.0, 0.0, 0.0, 0.0],  # DT
    [1.0, 0.0, 0.0, 0.0],  # A
    [1.0, 0.0, 0.0, 0.0],  # C
    [1.0, 0.0, 0.0, 0.0],  # G
    [1.0, 0.0, 0.0, 0.0],  # U
    [0.0, 0.0, 0.0, 0.0],  # UNK
    [1.0, 0.0, 0.0, 0.0],  # DN
    [1.0, 0.0, 0.0, 0.0],  # N
]

# The following chi angles are pi periodic: they can be rotated by a multiple
# of pi without affecting the structure.
chi_pi_periodic = [
    [0.0, 0.0, 0.0, 0.0],  # ALA
    [0.0, 0.0, 0.0, 0.0],  # ARG
    [0.0, 0.0, 0.0, 0.0],  # ASN
    [0.0, 1.0, 0.0, 0.0],  # ASP
    [0.0, 0.0, 0.0, 0.0],  # CYS
    [0.0, 0.0, 0.0, 0.0],  # GLN
    [0.0, 0.0, 1.0, 0.0],  # GLU
    [0.0, 0.0, 0.0, 0.0],  # GLY
    [0.0, 0.0, 0.0, 0.0],  # HIS
    [0.0, 0.0, 0.0, 0.0],  # ILE
    [0.0, 0.0, 0.0, 0.0],  # LEU
    [0.0, 0.0, 0.0, 0.0],  # LYS
    [0.0, 0.0, 0.0, 0.0],  # MET
    [0.0, 1.0, 0.0, 0.0],  # PHE
    [0.0, 0.0, 0.0, 0.0],  # PRO
    [0.0, 0.0, 0.0, 0.0],  # SER
    [0.0, 0.0, 0.0, 0.0],  # THR
    [0.0, 0.0, 0.0, 0.0],  # TRP
    [0.0, 1.0, 0.0, 0.0],  # TYR
    [0.0, 0.0, 0.0, 0.0],  # VAL
    [0.0, 0.0, 0.0, 0.0],  # DA
    [0.0, 0.0, 0.0, 0.0],  # DC
    [0.0, 0.0, 0.0, 0.0],  # DG
    [0.0, 0.0, 0.0, 0.0],  # DT
    [0.0, 0.0, 0.0, 0.0],  # A
    [0.0, 0.0, 0.0, 0.0],  # C
    [0.0, 0.0, 0.0, 0.0],  # G
    [0.0, 0.0, 0.0, 0.0],  # U
    [0.0, 0.0, 0.0, 0.0],  # UNK,
    [0.0, 0.0, 0.0, 0.0],  # DN,
    [0.0, 0.0, 0.0, 0.0],  # N,
]

# Proteins (AlphaFold2)
# Atoms positions relative to the 8 rigid groups, defined by the pre-omega, phi,
# psi and chi angles:
# 0: 'backbone group',
# 1: 'pre-omega-group', (empty)
# 2: 'phi-group', (currently empty, because it defines only hydrogens)
# 3: 'psi-group',
# 4,5,6,7: 'chi1,2,3,4-group'
# The atom positions are relative to the axis-end-atom of the corresponding
# rotation axis. The x-axis is in direction of the rotation axis, and the y-axis
# is defined such that the dihedral-angle-definiting atom (the last entry in
# chi_angles_atoms above) is in the xy-plane (with a positive y-coordinate).
# format: [atomname, group_idx, rel_position]

# Nucleotides
# 0: 'backbone group'
# 1: alpha
# 2: beta
# 3: gamma
# 4: delta
# 5: nu2
# 6: nu1
# 7: nu0
# 8: chi1

rigid_group_atom_positions = {
    "ALA": [
        ["N", 0, (-0.525, 1.363, 0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.526, -0.000, -0.000)],
        ["CB", 0, (-0.529, -0.774, -1.205)],
        ["O", 3, (0.627, 1.062, 0.000)],
    ],
    "ARG": [
        ["N", 0, (-0.524, 1.362, -0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.525, -0.000, -0.000)],
        ["CB", 0, (-0.524, -0.778, -1.209)],
        ["O", 3, (0.626, 1.062, 0.000)],
        ["CG", 4, (0.616, 1.390, -0.000)],
        ["CD", 5, (0.564, 1.414, 0.000)],
        ["NE", 6, (0.539, 1.357, -0.000)],
        ["NH1", 7, (0.206, 2.301, 0.000)],
        ["NH2", 7, (2.078, 0.978, -0.000)],
        ["CZ", 7, (0.758, 1.093, -0.000)],
    ],
    "ASN": [
        ["N", 0, (-0.536, 1.357, 0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.526, -0.000, -0.000)],
        ["CB", 0, (-0.531, -0.787, -1.200)],
        ["O", 3, (0.625, 1.062, 0.000)],
        ["CG", 4, (0.584, 1.399, 0.000)],
        ["ND2", 5, (0.593, -1.188, 0.001)],
        ["OD1", 5, (0.633, 1.059, 0.000)],
    ],
    "ASP": [
        ["N", 0, (-0.525, 1.362, -0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.527, 0.000, -0.000)],
        ["CB", 0, (-0.526, -0.778, -1.208)],
        ["O", 3, (0.626, 1.062, -0.000)],
        ["CG", 4, (0.593, 1.398, -0.000)],
        ["OD1", 5, (0.610, 1.091, 0.000)],
        ["OD2", 5, (0.592, -1.101, -0.003)],
    ],
    "CYS": [
        ["N", 0, (-0.522, 1.362, -0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.524, 0.000, 0.000)],
        ["CB", 0, (-0.519, -0.773, -1.212)],
        ["O", 3, (0.625, 1.062, -0.000)],
        ["SG", 4, (0.728, 1.653, 0.000)],
    ],
    "GLN": [
        ["N", 0, (-0.526, 1.361, -0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.526, 0.000, 0.000)],
        ["CB", 0, (-0.525, -0.779, -1.207)],
        ["O", 3, (0.626, 1.062, -0.000)],
        ["CG", 4, (0.615, 1.393, 0.000)],
        ["CD", 5, (0.587, 1.399, -0.000)],
        ["NE2", 6, (0.593, -1.189, -0.001)],
        ["OE1", 6, (0.634, 1.060, 0.000)],
    ],
    "GLU": [
        ["N", 0, (-0.528, 1.361, 0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.526, -0.000, -0.000)],
        ["CB", 0, (-0.526, -0.781, -1.207)],
        ["O", 3, (0.626, 1.062, 0.000)],
        ["CG", 4, (0.615, 1.392, 0.000)],
        ["CD", 5, (0.600, 1.397, 0.000)],
        ["OE1", 6, (0.607, 1.095, -0.000)],
        ["OE2", 6, (0.589, -1.104, -0.001)],
    ],
    "GLY": [
        ["N", 0, (-0.572, 1.337, 0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.517, -0.000, -0.000)],
        ["O", 3, (0.626, 1.062, -0.000)],
    ],
    "HIS": [
        ["N", 0, (-0.527, 1.360, 0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.525, 0.000, 0.000)],
        ["CB", 0, (-0.525, -0.778, -1.208)],
        ["O", 3, (0.625, 1.063, 0.000)],
        ["CG", 4, (0.600, 1.370, -0.000)],
        ["CD2", 5, (0.889, -1.021, 0.003)],
        ["ND1", 5, (0.744, 1.160, -0.000)],
        ["CE1", 5, (2.030, 0.851, 0.002)],
        ["NE2", 5, (2.145, -0.466, 0.004)],
    ],
    "ILE": [
        ["N", 0, (-0.493, 1.373, -0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.527, -0.000, -0.000)],
        ["CB", 0, (-0.536, -0.793, -1.213)],
        ["O", 3, (0.627, 1.062, -0.000)],
        ["CG1", 4, (0.534, 1.437, -0.000)],
        ["CG2", 4, (0.540, -0.785, -1.199)],
        ["CD1", 5, (0.619, 1.391, 0.000)],
    ],
    "LEU": [
        ["N", 0, (-0.520, 1.363, 0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.525, -0.000, -0.000)],
        ["CB", 0, (-0.522, -0.773, -1.214)],
        ["O", 3, (0.625, 1.063, -0.000)],
        ["CG", 4, (0.678, 1.371, 0.000)],
        ["CD1", 5, (0.530, 1.430, -0.000)],
        ["CD2", 5, (0.535, -0.774, 1.200)],
    ],
    "LYS": [
        ["N", 0, (-0.526, 1.362, -0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.526, 0.000, 0.000)],
        ["CB", 0, (-0.524, -0.778, -1.208)],
        ["O", 3, (0.626, 1.062, -0.000)],
        ["CG", 4, (0.619, 1.390, 0.000)],
        ["CD", 5, (0.559, 1.417, 0.000)],
        ["CE", 6, (0.560, 1.416, 0.000)],
        ["NZ", 7, (0.554, 1.387, 0.000)],
    ],
    "MET": [
        ["N", 0, (-0.521, 1.364, -0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.525, 0.000, 0.000)],
        ["CB", 0, (-0.523, -0.776, -1.210)],
        ["O", 3, (0.625, 1.062, -0.000)],
        ["CG", 4, (0.613, 1.391, -0.000)],
        ["SD", 5, (0.703, 1.695, 0.000)],
        ["CE", 6, (0.320, 1.786, -0.000)],
    ],
    "PHE": [
        ["N", 0, (-0.518, 1.363, 0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.524, 0.000, -0.000)],
        ["CB", 0, (-0.525, -0.776, -1.212)],
        ["O", 3, (0.626, 1.062, -0.000)],
        ["CG", 4, (0.607, 1.377, 0.000)],
        ["CD1", 5, (0.709, 1.195, -0.000)],
        ["CD2", 5, (0.706, -1.196, 0.000)],
        ["CE1", 5, (2.102, 1.198, -0.000)],
        ["CE2", 5, (2.098, -1.201, -0.000)],
        ["CZ", 5, (2.794, -0.003, -0.001)],
    ],
    "PRO": [
        ["N", 0, (-0.566, 1.351, -0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.527, -0.000, 0.000)],
        ["CB", 0, (-0.546, -0.611, -1.293)],
        ["O", 3, (0.621, 1.066, 0.000)],
        ["CG", 4, (0.382, 1.445, 0.0)],
        ["CD", 5, (0.477, 1.424, 0.0)],
    ],
    "SER": [
        ["N", 0, (-0.529, 1.360, -0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.525, -0.000, -0.000)],
        ["CB", 0, (-0.518, -0.777, -1.211)],
        ["O", 3, (0.626, 1.062, -0.000)],
        ["OG", 4, (0.503, 1.325, 0.000)],
    ],
    "THR": [
        ["N", 0, (-0.517, 1.364, 0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.526, 0.000, -0.000)],
        ["CB", 0, (-0.516, -0.793, -1.215)],
        ["O", 3, (0.626, 1.062, 0.000)],
        ["CG2", 4, (0.550, -0.718, -1.228)],
        ["OG1", 4, (0.472, 1.353, 0.000)],
    ],
    "TRP": [
        ["N", 0, (-0.521, 1.363, 0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.525, -0.000, 0.000)],
        ["CB", 0, (-0.523, -0.776, -1.212)],
        ["O", 3, (0.627, 1.062, 0.000)],
        ["CG", 4, (0.609, 1.370, -0.000)],
        ["CD1", 5, (0.824, 1.091, 0.000)],
        ["CD2", 5, (0.854, -1.148, -0.005)],
        ["CE2", 5, (2.186, -0.678, -0.007)],
        ["CE3", 5, (0.622, -2.530, -0.007)],
        ["NE1", 5, (2.140, 0.690, -0.004)],
        ["CH2", 5, (3.028, -2.890, -0.013)],
        ["CZ2", 5, (3.283, -1.543, -0.011)],
        ["CZ3", 5, (1.715, -3.389, -0.011)],
    ],
    "TYR": [
        ["N", 0, (-0.522, 1.362, 0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.524, -0.000, -0.000)],
        ["CB", 0, (-0.522, -0.776, -1.213)],
        ["O", 3, (0.627, 1.062, -0.000)],
        ["CG", 4, (0.607, 1.382, -0.000)],
        ["CD1", 5, (0.716, 1.195, -0.000)],
        ["CD2", 5, (0.713, -1.194, -0.001)],
        ["CE1", 5, (2.107, 1.200, -0.002)],
        ["CE2", 5, (2.104, -1.201, -0.003)],
        ["OH", 5, (4.168, -0.002, -0.005)],
        ["CZ", 5, (2.791, -0.001, -0.003)],
    ],
    "VAL": [
        ["N", 0, (-0.494, 1.373, -0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.527, -0.000, -0.000)],
        ["CB", 0, (-0.533, -0.795, -1.213)],
        ["O", 3, (0.627, 1.062, -0.000)],
        ["CG1", 4, (0.540, 1.429, -0.000)],
        ["CG2", 4, (0.533, -0.776, 1.203)],
    ],
    "DA": [
        ["OP1", 0, (-0.7319, 1.2920, 0.000)],
        ["P", 0, (0.000, 0.000, 0.000)],
        ["OP2", 0, (1.4855, 0.000, 0.000)],
        ["O5'", 0, (-0.4948, -0.8559, 1.2489)],
        ["C5'", 1, (0.7411, 1.2354, 0.000)],
        ["C4'", 2, (0.5207, 1.4178, 0.000)],
        ["C3'", 3, (0.6388, 1.3889, 0.000)],
        ["O4'", 3, (0.4804, -0.6610, -1.1947)],
        ["C1'", 6, (0.2998, 1.4906, 0.0000)],
        ["N9", 6, (1.4762, 1.8450, 0.7133)],
        ["C2'", 5, (0.4167, 1.4603, 0.0000)],
        ["O3'", 4, (0.4966, 1.3432, 0.0000)],
        ["C4", 8, (1.0591, 0.8359, 0.0000)],
        ["N3", 8, (1.1876, 1.9681, 0.6563)],
        ["C2", 8, (2.3468, 2.5190, 0.4196)],
        ["N1", 8, (3.3164, 2.0956, -0.3417)],
        ["C6", 8, (3.1559, 0.9546, -0.9869)],
        ["C5", 8, (1.9678, 0.2745, -0.8177)],
        ["N7", 8, (1.4954, -0.8947, -1.3248)],
        ["C8", 8, (0.3278, -1.0116, -0.8099)],
        ["N6", 8, (4.0967, 0.5227, -1.7225)],
    ],
    "DC": [
        ["OP1", 0, (-0.7319, 1.2920, 0.000)],
        ["P", 0, (0.000, 0.000, 0.000)],
        ["OP2", 0, (1.4855, 0.000, 0.000)],
        ["O5'", 0, (-0.4948, -0.8559, 1.2489)],
        ["C5'", 1, (0.7411, 1.2354, 0.000)],
        ["C4'", 2, (0.5207, 1.4178, 0.000)],
        ["C3'", 3, (0.6388, 1.3889, 0.000)],
        ["O4'", 3, (0.4804, -0.6610, -1.1947)],
        ["C1'", 6, (0.3291, 1.4833, 0.0000)],
        ["N1", 6, (1.5038, 1.8428, 0.7170)],
        ["C2'", 5, (0.4167, 1.4603, 0.0000)],
        ["O3'", 4, (0.4966, 1.3432, 0.0000)],
        ["C2", 8, (0.8923, 1.0332, 0.0000)],
        ["O2", 8, (0.6383, 2.0410, 0.6020)],
        ["N3", 8, (2.0166, 0.9037, -0.6631)],
        ["C4", 8, (2.2618, -0.1975, -1.3017)],
        ["N4", 8, (3.3803, -0.2807, -1.9426)],
        ["C5", 8, (1.3687, -1.2651, -1.3114)],
        ["C6", 8, (0.2615, -1.1249, -0.6529)],
    ],
    "DG": [
        ["OP1", 0, (-0.7319, 1.2920, 0.000)],
        ["P", 0, (0.000, 0.000, 0.000)],
        ["OP2", 0, (1.4855, 0.000, 0.000)],
        ["O5'", 0, (-0.4948, -0.8559, 1.2489)],
        ["C5'", 1, (0.7411, 1.2354, 0.000)],
        ["C4'", 2, (0.5207, 1.4178, 0.000)],
        ["C3'", 3, (0.6388, 1.3889, 0.000)],
        ["O4'", 3, (0.4804, -0.6610, -1.1947)],
        ["C1'", 6, (0.3226, 1.4925, 0.0000)],
        ["N9", 6, (1.4900, 1.8499, 0.6687)],
        ["C2'", 5, (0.4167, 1.4603, 0.0000)],
        ["O3'", 4, (0.4966, 1.3432, 0.0000)],
        ["C4", 8, (1.0086, 0.8849, 0.0000)],
        ["N3", 8, (1.0496, 2.0150, 0.6679)],
        ["C2", 8, (2.1496, 2.6645, 0.4679)],
        ["N1", 8, (3.1297, 2.2367, -0.3228)],
        ["C6", 8, (3.1039, 1.0771, -1.0158)],
        ["C5", 8, (1.9331, 0.3751, -0.8115)],
        ["N7", 8, (1.5201, -0.8107, -1.3144)],
        ["C8", 8, (0.3729, -0.9915, -0.8060)],
        ["N2", 8, (2.3472, 3.8096, 1.0562)],
        ["O6", 8, (4.0380, 0.7786, -1.7024)],
    ],
    "DT": [
        ["OP1", 0, (-0.7319, 1.2920, 0.000)],
        ["P", 0, (0.000, 0.000, 0.000)],
        ["OP2", 0, (1.4855, 0.000, 0.000)],
        ["O5'", 0, (-0.4948, -0.8559, 1.2489)],
        ["C5'", 1, (0.7411, 1.2354, 0.000)],
        ["C4'", 2, (0.5207, 1.4178, 0.000)],
        ["C3'", 3, (0.6388, 1.3889, 0.000)],
        ["O4'", 3, (0.4804, -0.6610, -1.1947)],
        ["C1'", 6, (0.3260, 1.4830, 0.0000)],
        ["N1", 6, (1.5172, 1.8567, 0.7052)],
        ["C2'", 5, (0.3435, 1.4761, 0.0000)],
        ["O3'", 4, (0.4966, 1.3432, 0.0000)],
        ["C2", 8, (0.9090, 0.9994, 0.0000)],
        ["O2", 8, (0.7438, 2.0419, 0.5454)],
        ["N3", 8, (2.0295, 0.7382, -0.6675)],
        ["C4", 8, (2.3265, -0.3977, -1.3208)],
        ["O4", 8, (3.3706, -0.5161, -1.8822)],
        ["C5", 8, (1.3322, -1.4021, -1.2827)],
        ["C7", 8, (1.5529, -2.6750, -1.9666)],
        ["C6", 8, (0.2324, -1.1610, -0.6324)],
    ],
    "A": [
        ["OP1", 0, (-0.7319, 1.2920, 0.000)],
        ["P", 0, (0.000, 0.000, 0.000)],
        ["OP2", 0, (1.4855, 0.000, 0.000)],
        ["O5'", 0, (-0.4948, -0.8559, 1.2489)],
        ["C5'", 1, (0.7289, 1.2185, 0.000)],
        ["C4'", 2, (0.5541, 1.4027, 0.000)],
        ["C3'", 3, (0.6673, 1.3669, 0.000)],
        ["O4'", 3, (0.4914, -0.6338, -1.2098)],
        ["C1'", 6, (0.2998, 1.4906, 0.0000)],
        ["N9", 6, (0.6843, 1.9988, 1.2090)],
        ["C2'", 5, (0.4641, 1.4573, 0.0000)],
        ["O2'", 7, (-1.9775, 1.3447, 0.0000)],
        ["O3'", 4, (0.5548, 1.3039, 0.000)],
        ["N1", 8, (2.1780, 3.2213, -0.2989)],
        ["C2", 8, (1.0879, 3.2302, 0.3745)],
        ["N3", 8, (0.2479, 2.2707, 0.5815)],
        ["C4", 8, (0.6207, 1.1816, 0.0000)],
        ["C5", 8, (1.7137, 1.0229, -0.7241)],
        ["C6", 8, (2.5253, 2.1104, -0.8708)],
        ["N6", 8, (3.6125, 2.1016, -1.5409)],
        ["N7", 8, (1.7832, -0.2340, -1.1740)],
        ["C8", 8, (0.7473, -0.7995, -0.7191)],
    ],
    "C": [
        ["OP1", 0, (-0.7319, 1.2920, 0.000)],
        ["P", 0, (0.000, 0.000, 0.000)],
        ["OP2", 0, (1.4855, 0.000, 0.000)],
        ["O5'", 0, (-0.4948, -0.8559, 1.2489)],
        ["C5'", 1, (0.7289, 1.2185, 0.000)],
        ["C4'", 2, (0.5541, 1.4027, 0.000)],
        ["C3'", 3, (0.6673, 1.3669, 0.000)],
        ["O4'", 3, (0.4914, -0.6338, -1.2098)],
        ["C1'", 6, (0.2993, 1.4958, 0.0000)],
        ["N1", 6, (0.6002, 2.0107, 1.2937)],
        ["C2'", 5, (0.4641, 1.4573, 0.0000)],
        ["O2'", 7, (-1.9703, 1.3476, 0.0000)],
        ["O3'", 4, (0.5548, 1.3039, 0.000)],
        ["C2", 8, (0.4599, 1.3075, -0.1217)],
        ["O2", 8, (-0.1636, 2.1569, 0.6203)],
        ["N3", 8, (1.5714, 1.6152, -0.6754)],
        ["C4", 8, (2.2175, 0.6759, -1.3322)],
        ["N4", 8, (3.3089, 1.0260, -1.9859)],
        ["C5", 8, (1.7732, -0.6640, -1.3484)],
        ["C6", 8, (0.6719, -0.9549, -0.6751)],
    ],
    "G": [
        ["OP1", 0, (-0.7319, 1.2920, 0.000)],
        ["P", 0, (0.000, 0.000, 0.000)],
        ["OP2", 0, (1.4855, 0.000, 0.000)],
        ["O5'", 0, (-0.4948, -0.8559, 1.2489)],
        ["C5'", 1, (0.7289, 1.2185, 0.000)],
        ["C4'", 2, (0.5541, 1.4027, 0.000)],
        ["C3'", 3, (0.6673, 1.3669, 0.000)],
        ["O4'", 3, (0.4914, -0.6338, -1.2098)],
        ["C1'", 6, (0.2976, 1.4937, 0.0000)],
        ["N9", 6, (0.5974, 1.9948, 1.2654)],
        ["C2'", 5, (0.4641, 1.4573, 0.0000)],
        ["O2'", 7, (-1.9759, 1.3450, 0.0000)],
        ["O3'", 4, (0.5548, 1.3039, 0.000)],
        ["N1", 8, (1.9303, 3.3536, -0.2684)],
        ["C2", 8, (0.8129, 3.3416, 0.3985)],
        ["N2", 8, (0.4814, 4.4795, 0.8986)],
        ["N3", 8, (0.0727, 2.2908, 0.5647)],
        ["C4", 8, (0.5463, 1.2335, 0.0000)],
        ["C5", 8, (1.6501, 1.1433, -0.6866)],
        ["C6", 8, (2.4354, 2.2766, -0.8568)],
        ["O6", 8, (3.4524, 2.3917, -1.4361)],
        ["N7", 8, (1.7981, -0.1202, -1.1107)],
        ["C8", 8, (0.7989, -0.7620, -0.6819)],
    ],
    "U": [
        ["OP1", 0, (-0.7319, 1.2920, 0.000)],
        ["P", 0, (0.000, 0.000, 0.000)],
        ["OP2", 0, (1.4855, 0.000, 0.000)],
        ["O5'", 0, (-0.4948, -0.8559, 1.2489)],
        ["C5'", 1, (0.7289, 1.2185, 0.000)],
        ["C4'", 2, (0.5541, 1.4027, 0.000)],
        ["C3'", 3, (0.6673, 1.3669, 0.000)],
        ["O4'", 3, (0.4914, -0.6338, -1.2098)],
        ["C1'", 6, (0.2976, 1.4688, 0.000)],
        ["N1", 6, (0.8326, 1.9765, 1.0916)],
        ["C2'", 5, (0.4641, 1.4573, 0.0000)],
        ["O2'", 7, (-1.9644, 1.3474, 0.0000)],
        ["O3'", 4, (0.5548, 1.3039, 0.0000)],
        ["C2", 8, (0.1738, 1.8286, 0.0000)],
        ["O2", 8, (-0.3051, 2.7218, 0.5543)],
        ["N3", 8, (1.2798, 1.9598, -0.6749)],
        ["C4", 8, (1.9070, 1.0235, -1.3364)],
        ["O4", 8, (2.8962, 1.2894, -1.8985)],
        ["C5", 8, (1.3118, -0.2053, -1.2930)],
        ["C6", 8, (0.2329, -0.3604, -0.6419)],
    ],
}
rigid_group_atom_positions["UNK"] = rigid_group_atom_positions["GLY"]
rigid_group_atom_positions["DN"] = rigid_group_atom_positions["DC"]
rigid_group_atom_positions["N"] = rigid_group_atom_positions["C"]
num_frames = 9

# atomc -> atom condensed
num_atomc = max([len(c) for c in restype3_to_atoms_index.values()])
restype_name_to_atomc_names = {}
for k in restype3_to_atoms:
    res_num_atoms = len(restype3_to_atoms)
    atom_names = restype3_to_atoms[k]
    if len(atom_names) < num_atomc:
        atom_names += [""] * (num_atomc - len(atom_names))
    restype_name_to_atomc_names[k] = atom_names


atomc_backbone_mask = np.zeros((1, num_atomc), dtype=np.float32)
atomc_backbone_mask[:, :4] = 1

atomc_names_arr = np.array(list(restype_name_to_atomc_names.values()), dtype=object)
element_names_arr = np.array(
    [
        [x if len(x) == 0 else x[:1] for x in y]
        for y in restype_name_to_atomc_names.values()
    ],
    dtype=object,
)


def _make_rigid_transformation_4x4(ex, ey, translation):
    """Create a rigid 4x4 transformation matrix from two axes and transl."""
    # Normalize ex.
    ex_normalized = ex / np.linalg.norm(ex)

    # make ey perpendicular to ex
    ey_normalized = ey - np.dot(ey, ex_normalized) * ex_normalized
    ey_normalized /= np.linalg.norm(ey_normalized)

    # compute ez as cross product
    eznorm = np.cross(ex_normalized, ey_normalized)
    m = np.stack([ex_normalized, ey_normalized, eznorm, translation]).transpose()
    m = np.concatenate([m, [[0.0, 0.0, 0.0, 1.0]]], axis=0)
    return m


# create an array with (restype, atomtype) --> rigid_group_idx
# and an array with (restype, atomtype, coord) for the atom positions
# and compute affine transformation matrices (4,4) from one rigid group to the
# previous group
restype_atomf_to_rigid_group = np.zeros(
    [full_num_residues, atom_type_num], dtype=np.int64
)
restype_atomf_mask = np.zeros([full_num_residues, atom_type_num], dtype=np.float32)
restype_atomf_rigid_group_positions = np.zeros(
    [full_num_residues, atom_type_num, 3], dtype=np.float32
)
restype_atomc_to_rigid_group = np.zeros(
    [full_num_residues, num_atomc], dtype=np.int64
)
restype_atomc_mask = np.zeros([full_num_residues, num_atomc], dtype=np.float32)
restype_atomc_rigid_group_positions = np.zeros(
    [full_num_residues, num_atomc, 3], dtype=np.float32
)
restype_atom3_rigid_group_positions = np.zeros(
    [full_num_residues, 3, 3], dtype=np.float32
)
restype_rigid_group_default_frame = np.zeros(
    [full_num_residues, num_frames, 4, 4], dtype=np.float32
)


def _make_rigid_group_constants():
    """Fill the arrays above."""
    for restype_letter in index_to_restype_1:
        restype = restype_1_to_index[restype_letter]
        resname = restype_1to3[restype_letter]
        for atomname, group_idx, atom_position in rigid_group_atom_positions[resname]:
            atomtype = atom_order[atomname]
            restype_atomf_to_rigid_group[restype, atomtype] = group_idx
            restype_atomf_mask[restype, atomtype] = 1
            restype_atomf_rigid_group_positions[restype, atomtype, :] = atom_position

            atomcidx = restype_name_to_atomc_names[resname].index(atomname)
            restype_atomc_to_rigid_group[restype, atomcidx] = group_idx
            restype_atomc_mask[restype, atomcidx] = 1
            restype_atomc_rigid_group_positions[restype, atomcidx, :] = atom_position

            if atomname in backbone_atoms:
                restype_atom3_rigid_group_positions[
                    restype, atomcidx, :
                ] = atom_position

        atom_positions = {
            name: np.array(pos) for name, _, pos in rigid_group_atom_positions[resname]
        }

        # Frame computations

        # backbone to backbone is the identity transform
        restype_rigid_group_default_frame[restype, 0, :, :] = np.eye(4)
        if restype3_is_prot(resname):
            # pre-omega-frame to backbone (currently dummy identity matrix)
            restype_rigid_group_default_frame[restype, 1, :, :] = np.eye(4)

            # phi-frame to backbone
            mat = _make_rigid_transformation_4x4(
                ex=atom_positions["N"] - atom_positions["CA"],
                ey=np.array([1.0, 0.0, 0.0]),
                translation=atom_positions["N"],
            )
            restype_rigid_group_default_frame[restype, 2, :, :] = mat

            # psi-frame to backbone
            mat = _make_rigid_transformation_4x4(
                ex=atom_positions["C"] - atom_positions["CA"],
                ey=atom_positions["CA"] - atom_positions["N"],
                translation=atom_positions["C"],
            )
            restype_rigid_group_default_frame[restype, 3, :, :] = mat

            # chi1-frame to backbone
            if chi_angles_mask[restype][0]:
                base_atom_names = chi_angles_atoms[resname][0]
                base_atom_positions = [atom_positions[name] for name in base_atom_names]
                mat = _make_rigid_transformation_4x4(
                    ex=base_atom_positions[2] - base_atom_positions[1],
                    ey=base_atom_positions[0] - base_atom_positions[1],
                    translation=base_atom_positions[2],
                )
                restype_rigid_group_default_frame[restype, 4, :, :] = mat

            # chi2-frame to chi1-frame
            # chi3-frame to chi2-frame
            # chi4-frame to chi3-frame
            # luckily all rotation axes for the next frame start at (0,0,0) of the
            # previous frame
            for chi_idx in range(1, 4):
                if chi_angles_mask[restype][chi_idx]:
                    axis_end_atom_name = chi_angles_atoms[resname][chi_idx][2]
                    axis_end_atom_position = atom_positions[axis_end_atom_name]
                    mat = _make_rigid_transformation_4x4(
                        ex=axis_end_atom_position,
                        ey=np.array([-1.0, 0.0, 0.0]),
                        translation=axis_end_atom_position,
                    )
                    restype_rigid_group_default_frame[restype, 4 + chi_idx, :, :] = mat

        elif restype3_is_na(resname):
            # alpha
            restype_rigid_group_default_frame[
                restype, 1, :, :
            ] = _make_rigid_transformation_4x4(
                ex=atom_positions["O5'"] - atom_positions["P"],
                ey=atom_positions["OP1"] - atom_positions["P"],
                translation=atom_positions["O5'"],
            )

            # beta
            restype_rigid_group_default_frame[
                restype, 2, :, :
            ] = _make_rigid_transformation_4x4(
                ex=atom_positions["C5'"],
                ey=np.array([-1.0, 0.0, 0.0]),
                translation=atom_positions["C5'"],
            )

            # gamma
            restype_rigid_group_default_frame[
                restype, 3, :, :
            ] = _make_rigid_transformation_4x4(
                ex=atom_positions["C4'"],
                ey=np.array([-1.0, 0.0, 0.0]),
                translation=atom_positions["C4'"],
            )

            # delta
            restype_rigid_group_default_frame[
                restype, 4, :, :
            ] = _make_rigid_transformation_4x4(
                ex=atom_positions["C3'"],
                ey=np.array([-1.0, 0.0, 0.0]),
                translation=atom_positions["C3'"],
            )

            # nu2
            restype_rigid_group_default_frame[
                restype, 5, :, :
            ] = _make_rigid_transformation_4x4(
                ex=atom_positions["C3'"],
                ey=np.array([-1.0, 0.0, 0.0]),
                translation=atom_positions["C3'"],
            )

            # nu1
            restype_rigid_group_default_frame[
                restype, 6, :, :
            ] = _make_rigid_transformation_4x4(
                ex=atom_positions["C2'"],
                ey=np.array([-1.0, 0.0, 0.0]),
                translation=atom_positions["C2'"],
            )

            # nu0
            restype_rigid_group_default_frame[
                restype, 7, :, :
            ] = _make_rigid_transformation_4x4(
                ex=atom_positions["C2'"],
                ey=np.array([-1.0, 0.0, 0.0]),
                translation=atom_positions["C2'"],
            )

            # chi1
            base_atom_names = chi_angles_atoms[resname][0]
            base_atom_positions = [atom_positions[name] for name in base_atom_names]
            mat = _make_rigid_transformation_4x4(
                ex=base_atom_positions[2],
                ey=np.array([-1.0, 0.0, 0.0]),
                translation=base_atom_positions[2],
            )
            restype_rigid_group_default_frame[restype, 8, :, :] = mat


_make_rigid_group_constants()

# Naming swaps for ambiguous atom names.
# Due to symmetries in the amino acids the naming of atoms is ambiguous in
# 4 of the 20 amino acids.
# (The LDDT paper lists 7 amino acids as ambiguous, but the naming ambiguities
# in LEU, VAL and ARG can be resolved by using the 3d constellations of
# the 'ambiguous' atoms and their neighbours)
residue_atom_renaming_swaps = {
    "ASP": {"OD1": "OD2"},
    "GLU": {"OE1": "OE2"},
    "PHE": {"CD1": "CD2", "CE1": "CE2"},
    "TYR": {"CD1": "CD2", "CE1": "CE2"},
}

# Van der Waals radii [Angstrom] of the atoms (from Wikipedia)
van_der_waals_radius = {
    "C": 1.7,
    "N": 1.55,
    "O": 1.52,
    "S": 1.8,
    "P": 1.8,
    "": 0.0,
}


def sequence_to_onehot(
    sequence: str, mapping, map_unknown_to_x: bool = False
) -> torch.LongTensor:
    """Maps the given sequence into a one-hot encoded matrix.
    Args:
      sequence: An amino acid sequence.
      mapping: A dictionary mapping amino acids to integers.
      map_unknown_to_x: If True, any amino acid that is not in the mapping will be
        mapped to the unknown amino acid 'X'. If the mapping doesn't contain
        amino acid 'X', an error will be thrown. If False, any amino acid not in
        the mapping will throw an error.
    Returns:
      A numpy array of shape (seq_len, num_unique_aas) with one-hot encoding of
      the sequence.
    Raises:
      ValueError: If the mapping doesn't contain values from 0 to
        num_unique_aas - 1 without any gaps.
    """
    num_entries = max(mapping.values()) + 1

    if sorted(set(mapping.values())) != list(range(num_entries)):
        raise ValueError(
            "The mapping must have values from 0 to num_unique_aas-1 "
            "without any gaps. Got: %s" % sorted(mapping.values())
        )

    one_hot_arr = torch.zeros(len(sequence), num_entries, dtype=torch.long)

    for aa_index, aa_type in enumerate(sequence):
        if map_unknown_to_x:
            if aa_type.isalpha() and aa_type.isupper():
                aa_id = mapping.get(aa_type, mapping["X"])
            else:
                raise ValueError(f"Invalid character in the sequence: {aa_type}")
        else:
            aa_id = mapping[aa_type]
        one_hot_arr[aa_index, aa_id] = 1

    return one_hot_arr


def get_chi_atom_indices():
    """Returns atom indices needed to compute chi angles for all residue types.
    Returns:
      A tensor of shape [residue_types=21, chis=4, atoms=4]. The residue types are
      in the order specified in residue_constants.restypes + unknown residue type
      at the end. For chi angles which are not defined on the residue, the
      positions indices are by default set to 0.
    """
    chi_atom_indices = []
    for residue_name in index_to_restype_1:
        residue_name = restype_1to3[residue_name]
        residue_chi_angles = chi_angles_atoms[residue_name]
        atom_indices = []
        for chi_angle in residue_chi_angles:
            atom_indices.append([atom_order[atom] for atom in chi_angle])
        for _ in range(4 - len(atom_indices)):
            atom_indices.append([0, 0, 0, 0])  # For chi angles not defined on the AA.
        chi_atom_indices.append(atom_indices)

    chi_atom_indices.append([[0, 0, 0, 0]] * 4)  # For UNKNOWN residue.

    return np.asarray(chi_atom_indices)


chi_atom_indices = get_chi_atom_indices()

Bond = namedtuple("Bond", ["atom1_name", "atom2_name", "length", "stddev"])
BondAngle = namedtuple(
    "BondAngle", ["atom1_name", "atom2_name", "atom3name", "angle_rad", "stddev"]
)


@functools.lru_cache(maxsize=None)
def load_stereo_chemical_props() -> Tuple[
    Mapping[str, List[Bond]], Mapping[str, List[Bond]], Mapping[str, List[BondAngle]]
]:
    """Load stereo_chemical_props.txt into a nice structure.
    Load literature values for bond lengths and bond angles and translate
    bond angles into the length of the opposite edge of the triangle
    ("residue_virtual_bonds").
    Returns:
      residue_bonds: Dict that maps resname -> list of Bond tuples.
      residue_virtual_bonds: Dict that maps resname -> list of Bond tuples.
      residue_bond_angles: Dict that maps resname -> list of BondAngle tuples.
    """
    stereo_chemical_props_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "stereo_chemical_props.txt"
    )
    with open(stereo_chemical_props_path, "rt") as f:
        stereo_chemical_props = f.read()
    lines_iter = iter(stereo_chemical_props.splitlines())
    # Load bond lengths.
    residue_bonds = {}
    next(lines_iter)  # Skip header line.
    for line in lines_iter:
        if line.strip() == "-":
            break
        bond, resname, length, stddev = line.split()
        atom1, atom2 = bond.split("-")
        if resname not in residue_bonds:
            residue_bonds[resname] = []
        residue_bonds[resname].append(Bond(atom1, atom2, float(length), float(stddev)))
    residue_bonds["UNK"] = []

    # Load bond angles.
    residue_bond_angles = {}
    next(lines_iter)  # Skip empty line.
    next(lines_iter)  # Skip header line.
    for line in lines_iter:
        if line.strip() == "-":
            break
        bond, resname, angle_degree, stddev_degree = line.split()
        atom1, atom2, atom3 = bond.split("-")
        if resname not in residue_bond_angles:
            residue_bond_angles[resname] = []
        residue_bond_angles[resname].append(
            BondAngle(
                atom1,
                atom2,
                atom3,
                float(angle_degree) / 180.0 * np.pi,
                float(stddev_degree) / 180.0 * np.pi,
            )
        )
    residue_bond_angles["UNK"] = []

    def make_bond_key(atom1_name, atom2_name):
        """Unique key to lookup bonds."""
        return "-".join(sorted([atom1_name, atom2_name]))

    # Translate bond angles into distances ("virtual bonds").
    residue_virtual_bonds = {}
    for resname, bond_angles in residue_bond_angles.items():
        # Create a fast lookup dict for bond lengths.
        bond_cache = {}
        for b in residue_bonds[resname]:
            bond_cache[make_bond_key(b.atom1_name, b.atom2_name)] = b
        residue_virtual_bonds[resname] = []
        for ba in bond_angles:
            bond1 = bond_cache[make_bond_key(ba.atom1_name, ba.atom2_name)]
            bond2 = bond_cache[make_bond_key(ba.atom2_name, ba.atom3name)]

            # Compute distance between atom1 and atom3 using the law of cosines
            # c^2 = a^2 + b^2 - 2ab*cos(gamma).
            gamma = ba.angle_rad
            length = np.sqrt(
                bond1.length ** 2
                + bond2.length ** 2
                - 2 * bond1.length * bond2.length * np.cos(gamma)
            )

            # Propagation of uncertainty assuming uncorrelated errors.
            dl_outer = 0.5 / length
            dl_dgamma = (2 * bond1.length * bond2.length * np.sin(gamma)) * dl_outer
            dl_db1 = (2 * bond1.length - 2 * bond2.length * np.cos(gamma)) * dl_outer
            dl_db2 = (2 * bond2.length - 2 * bond1.length * np.cos(gamma)) * dl_outer
            stddev = np.sqrt(
                (dl_dgamma * ba.stddev) ** 2
                + (dl_db1 * bond1.stddev) ** 2
                + (dl_db2 * bond2.stddev) ** 2
            )
            residue_virtual_bonds[resname].append(
                Bond(ba.atom1_name, ba.atom3name, length, stddev)
            )

    return (residue_bonds, residue_virtual_bonds, residue_bond_angles)


def get_atomc_dists_bounds(overlap_tolerance=1.5, bond_length_tolerance_factor=15):
    """compute upper and lower bounds for bonds to assess violations."""
    restype_atomc_bond_lower_bound = np.zeros(
        [full_num_residues, num_atomc, num_atomc], np.float32
    )
    restype_atomc_bond_upper_bound = np.zeros(
        [full_num_residues, num_atomc, num_atomc], np.float32
    )
    restype_atomc_bond_stddev = np.zeros(
        [full_num_residues, num_atomc, num_atomc], np.float32
    )
    residue_bonds, residue_virtual_bonds, _ = load_stereo_chemical_props()
    for restype, restype_letter in enumerate(index_to_restype_1):
        resname = restype_1to3[restype_letter]
        atom_list = restype_name_to_atomc_names[resname]

        # create lower and upper bounds for clashes
        for atom1_idx, atom1_name in enumerate(atom_list):
            if not atom1_name:
                continue
            atom1_radius = van_der_waals_radius[atom1_name[0]]
            for atom2_idx, atom2_name in enumerate(atom_list):
                if (not atom2_name) or atom1_idx == atom2_idx:
                    continue
                atom2_radius = van_der_waals_radius[atom2_name[0]]
                lower = atom1_radius + atom2_radius - overlap_tolerance
                upper = 1e10
                restype_atomc_bond_lower_bound[restype, atom1_idx, atom2_idx] = lower
                restype_atomc_bond_lower_bound[restype, atom2_idx, atom1_idx] = lower
                restype_atomc_bond_upper_bound[restype, atom1_idx, atom2_idx] = upper
                restype_atomc_bond_upper_bound[restype, atom2_idx, atom1_idx] = upper

        if resname in residue_bonds:
            # overwrite lower and upper bounds for bonds and angles
            for b in residue_bonds[resname] + residue_virtual_bonds[resname]:
                atom1_idx = atom_list.index(b.atom1_name)
                atom2_idx = atom_list.index(b.atom2_name)
                lower = b.length - bond_length_tolerance_factor * b.stddev
                upper = b.length + bond_length_tolerance_factor * b.stddev
                restype_atomc_bond_lower_bound[restype, atom1_idx, atom2_idx] = lower
                restype_atomc_bond_lower_bound[restype, atom2_idx, atom1_idx] = lower
                restype_atomc_bond_upper_bound[restype, atom1_idx, atom2_idx] = upper
                restype_atomc_bond_upper_bound[restype, atom2_idx, atom1_idx] = upper
                restype_atomc_bond_stddev[restype, atom1_idx, atom2_idx] = b.stddev
                restype_atomc_bond_stddev[restype, atom2_idx, atom1_idx] = b.stddev
    return {
        "lower_bound": restype_atomc_bond_lower_bound,  # shape (28,23,23)
        "upper_bound": restype_atomc_bond_upper_bound,  # shape (28,23,23)
        "stddev": restype_atomc_bond_stddev,  # shape (28,23,23)
    }


atomc_dists_bounds = get_atomc_dists_bounds()

# Between-residue bond lengths for general bonds (first element) and for Proline
# (second element).
between_res_bond_length_c_n = [1.329, 1.341]
between_res_bond_length_stddev_c_n = [0.014, 0.016]

# Between-residue cos_angles.
between_res_cos_angles_c_n_ca = [-0.5203, 0.0353]  # degrees: 121.352 +- 2.315
between_res_cos_angles_ca_c_n = [-0.4473, 0.0311]  # degrees: 116.568 +- 1.995


def rename_aatype_to_convention(aatype):
    if torch.is_tensor(aatype):
        new_aatype = aatype.clone()
    else:
        new_aatype = aatype.copy()
    new_aatype[new_aatype == restype_3_to_index["UNK"]] = restype_3_to_index["GLY"]
    new_aatype[new_aatype == restype_3_to_index["DN"]] = restype_3_to_index["DC"]
    new_aatype[new_aatype == restype_3_to_index["N"]] = restype_3_to_index["C"]
    return new_aatype


def select_torsion_angles(input, aatype):
    new_aatype = rename_aatype_to_convention(aatype)
    chi_angles = einops.rearrange(
        input[..., 3:, :],
        "... (f a) d -> ... f d a",
        f=5,
        a=canonical_num_residues,
        d=2,
    )[torch.arange(len(new_aatype)), ..., new_aatype]
    input_torsion_angles = torch.cat((input[..., :3, :], chi_angles), dim=-2)
    input_torsion_angles = torch.nn.functional.normalize(input_torsion_angles, dim=-1)
    return input_torsion_angles


num_net_torsions = canonical_num_residues * 5 + 3

nuc_torsion_frames = [
    ["OP1", "P", "O5'", "C5'"],  # alpha, 1
    ["P", "O5'", "C5'", "C4'"],  # beta, 2
    ["O5'", "C5'", "C4'", "C3'"],  # gamma, 3
    ["C5'", "C4'", "C3'", "O3'"],  # delta, 4
    ["C5'", "C4'", "C3'", "C2'"],  # nu2, 5
    ["C4'", "C3'", "C2'", "C1'"],  # nu1, 6
    ["C3'", "C2'", "C1'", "O2'"],  # nu0, 7
]


nuc_torsion_atom_indices = []
nuc_torsion_atom_mask = []
for resname in ["DA", "DC", "DG", "DT", "A", "C", "G", "U", "DN", "N"]:
    resname_torsion_atom_indices = []
    resname_torsion_atom_mask = []
    for torsion_atoms in nuc_torsion_frames[:-1]:
        resname_torsion_atom_indices.append(
            [atom_order[atom] for atom in torsion_atoms]
        )
        resname_torsion_atom_mask.append(1)
    if len(resname) == 2:
        # Is DNA
        resname_torsion_atom_indices.append([0, 0, 0, 0])
        resname_torsion_atom_mask.append(0)
    else:
        resname_torsion_atom_indices.append(
            [atom_order[atom] for atom in nuc_torsion_frames[-1]]
        )
        resname_torsion_atom_mask.append(1)
    nuc_torsion_atom_indices.append(resname_torsion_atom_indices)
    nuc_torsion_atom_mask.append(resname_torsion_atom_mask)
nuc_torsion_atom_indices = np.array(nuc_torsion_atom_indices, dtype=np.int64)
nuc_torsion_atom_mask = np.array(nuc_torsion_atom_mask, dtype=np.int64)
