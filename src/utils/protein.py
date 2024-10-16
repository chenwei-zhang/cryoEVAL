######################################################
### Credit from ModelAngelo evaluation functions  ####
######################################################


import dataclasses
import pickle
import warnings
from typing import Dict, List

import numpy as np
import torch
from Bio.PDB import MMCIFParser, PDBParser

import utils.residue_constants as _rc
from utils.affine_utils import (
    affine_composition,
    affine_from_3_points,
    affine_from_tensor4x4,
    affine_mul_rots,
    affine_mul_vecs,
    fill_rotation_matrix,
    invert_affine,
)

# Complete sequence of chain IDs supported by the PDB format.
PDB_CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
PDB_MAX_CHAINS = len(PDB_CHAIN_IDS)  # := 62.

PROTEIN_KEYS = [
    "atom_positions",
    "atomc_positions",
    "aatype",
    "atom_mask",
    "atomc_mask",
    "residue_index",
    "chain_index",
    "chain_id",
    "b_factors",
    "rigidgroups_gt_frames",
    "rigidgroups_gt_exists",
    "rigidgroups_group_exists",
    "rigidgroups_group_is_ambiguous",
    "rigidgroups_alt_gt_frames",
    "torsion_angles_sin_cos",
    "alt_torsion_angles_sin_cos",
    "torsion_angles_mask",
    "unified_seq",
    "unified_seq_len",
    "residue_to_seq_id",
    "residue_to_lm_embedding",
    "chain_idx_to_residues",
    "prot_mask",
]


@dataclasses.dataclass(frozen=False)
class Protein:
    """Protein structure representation."""

    # Cartesian coordinates of atoms in angstroms. The atom types correspond to
    # _rc.atom_types, i.e. the first three are N, CA, CB.
    atom_positions: np.ndarray  # [num_res, 65, 3]

    # Cartesian coordinates of atoms in angstroms. The atom types correspond to
    # _rc.atom_types, i.e. the first three are N, CA, CB.
    atomc_positions: np.ndarray  # [num_res, _rc.num_atomc, 3]

    # Amino-acid type for each residue represented as an integer between 0 and
    # 20, where 20 is 'X'.
    aatype: np.ndarray  # [num_res]

    # Binary float mask to indicate presence of a particular atom. 1.0 if an atom
    # is present and 0.0 if not. This should be used for loss masking.
    atom_mask: np.ndarray  # [num_res, num_atom_type]

    # Same as above, but for atomc
    atomc_mask: np.ndarray

    # Residue index as used in PDB. It is not necessarily continuous or 0-indexed.
    residue_index: np.ndarray  # [num_res]

    # 0-indexed number corresponding to the chain in the protein that this residue
    # belongs to.
    chain_index: np.ndarray  # [num_res]

    # The original Chain ID string list that the chain_indices correspond to.
    chain_id: np.ndarray  # [num_chains]

    # B-factors, or temperature factors, of each residue (in sq. angstroms units),
    # representing the displacement of the residue from its ground truth mean
    # value.
    b_factors: np.ndarray  # [num_res, num_atom_type]

    # 8 Frames corresponding to 'atom_positions'
    # represented as affines.
    rigidgroups_gt_frames: np.ndarray  # (num_res, 8, 3, 4)

    # Mask denoting whether the atom positions for
    # the given frame are available in the ground truth, e.g. if they were
    # resolved in the experiment.
    rigidgroups_gt_exists: np.ndarray  # (num_res, 8)

    # Mask denoting whether given group is in
    # principle present for given amino acid type.
    rigidgroups_group_exists: np.ndarray  # (num_res, 8)

    # Mask denoting whether frame is
    # affected by naming ambiguity.
    rigidgroups_group_is_ambiguous: np.ndarray  # (num_res, 8)

    # 8 Frames with alternative atom renaming
    # corresponding to 'all_atom_positions' represented as affines
    rigidgroups_alt_gt_frames: np.ndarray  # (num_res, 8, 3, 4)

    # Array where the final
    # 2 dimensions denote sin and cos respectively
    torsion_angles_sin_cos: np.ndarray  # (num_res, 7, 2)

    # Same as 'torsion_angles_sin_cos', but
    # with the angle shifted by pi for all chi angles affected by the naming
    # ambiguities.
    alt_torsion_angles_sin_cos: np.ndarray  # (num_res, 7, 2)

    # Mask for which chi angles are present.
    torsion_angles_mask: np.ndarray  # (num_res, 7)

    # Unique sequences concatenated together, the ||| token is the delimiter
    unified_seq: str

    # Length of unified_seq, different than len(unified_seq) because of ||| tokens
    unified_seq_len: int

    # The mapping of residues to where they are in the unified_seq
    residue_to_seq_id: np.ndarray  # (num_res,)

    # The mapping of residues to their language model static embeddings, can be empty at construction
    residue_to_lm_embedding: np.ndarray  # (unified_seq_len, embedding_dim)

    # Chains to residues
    chain_idx_to_residues: List[np.ndarray]
    # Whether or not the residue is a protein residue
    prot_mask: np.ndarray  # (num_res,)

    keys = PROTEIN_KEYS


def get_protein_from_file_path(file_path: str, chain_id: str = None) -> Protein:
    """Takes a file path containing a PDB/mmCIF file and constructs a Protein object.
    WARNING: All non-standard residue types will be ignored. All
      non-standard atoms will be ignored.
    Args:
      pdb_str: The path to the PDB file
      chain_id: If chain_id is specified (e.g. A), then only that chain
        is parsed. Otherwise all chains are parsed.
    Returns:
      A new `Protein` parsed from the pdb contents.
    """
    if file_path.split(".")[-1][-3:] == "pdb":
        parser = PDBParser(QUIET=True)
    elif file_path.split(".")[-1][-3:] == "cif":
        parser = MMCIFParser(QUIET=True)
    else:
        raise RuntimeError("Unknown type for structure file:", file_path[-3:])
    structure = parser.get_structure("none", file_path)
    models = list(structure.get_models())
    model = models[0]

    atom_positions = []
    atomc_positions = []
    aatype = []
    atom_mask = []
    atomc_mask = []
    residue_index = []
    chain_ids = []
    chain_idx_to_residues = []
    b_factors = []

    # Sequence related
    unified_seq = []
    residue_to_seq_id = []
    temp_sequences_seen = {}
    seq_len_so_far = 0
    residue_count = 0

    for chain in model:
        if chain_id is not None and chain.id != chain_id:
            continue
        chain_seq = []
        chain_res_ids = []
        chain_aatype = []
        for res in chain:
            if res.resname not in _rc.restype_3to1:
                continue
            if res.id[2] != " ":
                warnings.warn(
                    f"PDB contains an insertion code at chain {chain.id} and residue "
                    f"index {res.id[1]}, {res.id[2]}. These are not supported."
                )
            res_shortname = _rc.restype_3to1[res.resname]
            restype_idx = _rc.restype_order.get(res_shortname, _rc.restype_num)

            pos = np.zeros((_rc.atom_type_num, 3))
            posc = np.zeros((_rc.num_atomc, 3))
            mask = np.zeros((_rc.atom_type_num,))
            maskc = np.zeros((_rc.num_atomc,))
            res_b_factors = np.zeros((_rc.atom_type_num,))
            for atom in res:
                if atom.name not in _rc.restype3_to_atoms_index[res.resname]:
                    continue
                pos[_rc.atom_order[atom.name]] = atom.coord
                posc[_rc.restype3_to_atoms_index[res.resname][atom.name]] = atom.coord
                mask[_rc.atom_order[atom.name]] = 1.0
                maskc[_rc.restype3_to_atoms_index[res.resname][atom.name]] = 1.0
                res_b_factors[_rc.atom_order[atom.name]] = atom.bfactor
            if np.sum(mask) < 0.5:
                # If no known atom positions are reported for the residue then skip it.
                continue
            if res_shortname.isupper():
                chain_seq.append(res_shortname)
            chain_aatype.append(restype_idx)
            atom_positions.append(pos)
            atomc_positions.append(posc)
            atom_mask.append(mask)
            atomc_mask.append(maskc)
            residue_index.append(res.id[1])
            chain_ids.append(chain.id)
            b_factors.append(res_b_factors)
            chain_res_ids.append(residue_count)
            residue_count += 1
        aatype.extend(chain_aatype)
        chain_idx_to_residues.append(np.array(chain_res_ids, dtype=np.int64))
        chain_aatype = np.array(chain_aatype, dtype=int)
        if len(chain_seq) > 0:
            chain_seq = "".join(chain_seq)
            if chain_seq not in temp_sequences_seen:
                temp_sequences_seen[chain_seq] = seq_len_so_far
                chain_residue_to_seq_id = list(
                    range(seq_len_so_far, seq_len_so_far + len(chain_seq))
                )
                unified_seq.append(chain_seq)
                seq_len_so_far += len(chain_seq)
            else:
                offset = temp_sequences_seen[chain_seq]
                chain_residue_to_seq_id = list(range(offset, offset + len(chain_seq)))
        else:
            chain_residue_to_seq_id = []
        if len(chain_seq) < len(chain_aatype):
            res_to_seq_id_array = np.array(chain_residue_to_seq_id, dtype=int)
            chain_residue_to_seq_id = np.zeros_like(chain_aatype)
            chain_residue_to_seq_id[chain_aatype < _rc.num_prot] = res_to_seq_id_array
            chain_residue_to_seq_id[chain_aatype >= _rc.num_prot] = -1
        residue_to_seq_id.extend(chain_residue_to_seq_id)
    if len(unified_seq) > 0:
        unified_seq = "|||".join(unified_seq) if len(unified_seq) > 1 else unified_seq[0]
    unified_seq_len = seq_len_so_far
    residue_to_seq_id = np.array(residue_to_seq_id, dtype=int)
    # Chain IDs are usually characters so map these to ints.
    unique_chain_ids = np.unique(chain_ids)
    chain_id_mapping = {cid: n for n, cid in enumerate(unique_chain_ids)}
    chain_index = np.array([chain_id_mapping[cid] for cid in chain_ids])

    atom_positions = np.array(atom_positions)
    atomc_positions = np.array(atomc_positions)
    atom_mask = np.array(atom_mask)
    atomc_mask = np.array(atomc_mask)
    aatype = np.array(aatype)
    residue_index = np.array(residue_index)
    b_factors = np.array(b_factors)
    frames = atomf_to_frames(
        aatype=aatype, all_atom_positions=atom_positions, all_atom_mask=atom_mask
    )
    torsion_angles = atomf_to_torsion_angles(
        aatype=aatype[None],
        all_atom_positions=atom_positions[None],
        all_atom_mask=atom_mask[None],
    )

    return Protein(
        atom_positions=atom_positions,
        atomc_positions=atomc_positions,
        atom_mask=atom_mask,
        atomc_mask=atomc_mask,
        aatype=aatype,
        residue_index=residue_index,
        chain_index=chain_index,
        chain_id=unique_chain_ids,
        b_factors=b_factors,
        unified_seq=unified_seq,
        unified_seq_len=unified_seq_len,
        residue_to_seq_id=residue_to_seq_id,
        residue_to_lm_embedding=None,
        chain_idx_to_residues=chain_idx_to_residues,
        prot_mask=aatype < _rc.num_prot,
        **frames,
        **torsion_angles,
    )


def get_protein_empty_except(**kwargs) -> Protein:
    protein_dict = dict.fromkeys(PROTEIN_KEYS)
    for (k, v) in kwargs.items():
        protein_dict[k] = v
    return Protein(**protein_dict)


def atomf_to_frames(
    aatype: np.ndarray,  # (...)
    all_atom_positions: np.ndarray,  # (..., 37, 3)
    all_atom_mask: np.ndarray,  # (..., 37)
) -> Dict[str, np.ndarray]:
    """Computes the frames for the up to 9 rigid groups for each residue.
    The rigid groups are defined by the possible torsions in a given amino acid.
    We group the atoms according to their dependence on the torsion angles into
    "rigid groups".  E.g., the position of atoms in the chi2-group depend on
    chi1 and chi2, but do not depend on chi3 or chi4.
    Jumper et al. (2021) Suppl. Table 2 and corresponding text.
    Args:
      aatype: Amino acid type, given as array with integers.
      all_atom_positions: atomf representation of all atom coordinates.
      all_atom_mask: atomf representation of mask on all atom coordinates.
    Returns:
      Dictionary containing:
        * 'rigidgroups_gt_frames': 9 Frames corresponding to 'all_atom_positions'
             represented as flat 12 dimensional array.
        * 'rigidgroups_gt_exists': Mask denoting whether the atom positions for
            the given frame are available in the ground truth, e.g. if they were
            resolved in the experiment.
        * 'rigidgroups_group_exists': Mask denoting whether given group is in
            principle present for given amino acid type.
        * 'rigidgroups_group_is_ambiguous': Mask denoting whether frame is
            affected by naming ambiguity.
        * 'rigidgroups_alt_gt_frames': 9 Frames with alternative atom renaming
            corresponding to 'all_atom_positions' represented as flat
            12 dimensional array.
    """
    # Proteins:
    # 0: 'backbone group',
    # 1: 'pre-omega-group', (empty)
    # 2: 'phi-group', (currently empty, because it defines only hydrogens)
    # 3: 'psi-group',
    # 4,5,6,7: 'chi1,2,3,4-group'
    aatype_in_shape = aatype.shape

    # If there is a batch axis, just flatten it away, and reshape everything
    # back at the end of the function.
    aatype = np.reshape(aatype, [-1])
    all_atom_positions = np.reshape(all_atom_positions, [-1, _rc.num_atoms, 3])
    all_atom_mask = np.reshape(all_atom_mask, [-1, _rc.num_atoms])
    N = len(aatype)

    # Create an array with the atom names.
    # shape (num_restypes, num_rigidgroups, 3_atoms): (28, 9, 3)
    restype_rigidgroup_base_atom_names = np.full(
        [_rc.full_num_residues, _rc.num_frames, 3], "", dtype=object
    )

    # 4,5,6,7: 'chi1,2,3,4-group'
    for restype, restype_letter in enumerate(_rc.index_to_restype_1):
        resname = _rc.restype_1to3[restype_letter]
        if _rc.restype3_is_prot(resname):
            # 0: backbone frame
            restype_rigidgroup_base_atom_names[restype, 0, :] = ["C", "CA", "N"]
            # 3: 'psi-group'
            restype_rigidgroup_base_atom_names[restype, 3, :] = ["CA", "C", "O"]
            for chi_idx in range(4):
                if _rc.chi_angles_mask[restype][chi_idx]:
                    atom_names = _rc.chi_angles_atoms[resname][chi_idx]
                    restype_rigidgroup_base_atom_names[
                        restype, chi_idx + 4, :
                    ] = atom_names[1:]
        else:
            # 0: backbone frame
            restype_rigidgroup_base_atom_names[restype, 0, :] = ["OP2", "P", "OP1"]
            # 1: alpha
            restype_rigidgroup_base_atom_names[restype, 1, :] = ["P", "O5'", "C5'"]
            # 2: beta
            restype_rigidgroup_base_atom_names[restype, 2, :] = ["O5'", "C5'", "C4'"]
            # 3: gamma
            restype_rigidgroup_base_atom_names[restype, 3, :] = ["C5'", "C4'", "C3'"]
            # 4: delta
            restype_rigidgroup_base_atom_names[restype, 4, :] = ["C4'", "C3'", "O3'"]
            # 5: nu2
            restype_rigidgroup_base_atom_names[restype, 5, :] = ["C4'", "C3'", "C2'"]
            # 6: nu1
            restype_rigidgroup_base_atom_names[restype, 6, :] = ["C3'", "C2'", "C1'"]
            # 7: nu0
            restype_rigidgroup_base_atom_names[restype, 7, :] = ["C2'", "C1'", "O2'"]
            # 8: chi1
            restype_rigidgroup_base_atom_names[restype, 8, :] = _rc.chi_angles_atoms[
                resname
            ][0][1:]

    # Create mask for existing rigid groups.
    restype_rigidgroup_mask = np.zeros(
        [_rc.full_num_residues, _rc.num_frames], dtype=np.float32
    )
    restype_rigidgroup_mask[:, 0] = 1
    restype_rigidgroup_mask[:, 3] = 1
    restype_rigidgroup_mask[: _rc.num_prot, 4:-1] = _rc.chi_angles_mask[: _rc.num_prot]
    restype_rigidgroup_mask[_rc.num_prot :, :] = 1
    restype_rigidgroup_mask[
        _rc.num_prot : _rc.num_prot + 4, 7
    ] = 0  # DNA does not have nu0

    # Translate atom names into atomf indices.
    lookuptable = _rc.atom_order.copy()
    lookuptable[""] = 0
    restype_rigidgroup_base_atomf_idx = np.vectorize(lambda x: lookuptable[x])(
        restype_rigidgroup_base_atom_names
    )

    # Compute the gather indices for all residues in the chain.
    # shape (N, 9, 3)
    residx_rigidgroup_base_atomf_idx = restype_rigidgroup_base_atomf_idx[aatype]

    # Gather the base atom positions for each rigid group.
    # Resulting shape: N, 9, 3, 3
    base_atom_pos_idx = (
        residx_rigidgroup_base_atomf_idx
        + np.arange(N * _rc.num_atoms, step=_rc.num_atoms)[..., None, None]
    )
    base_atom_pos = np.take(
        all_atom_positions.reshape(-1, 3), base_atom_pos_idx, axis=0
    )

    # Compute the Affines.
    gt_frames = affine_from_3_points(
        point_on_neg_x_axis=torch.Tensor(base_atom_pos[:, :, 0, :]),
        origin=torch.Tensor(base_atom_pos[:, :, 1, :]),
        point_on_xy_plane=torch.Tensor(base_atom_pos[:, :, 2, :]),
    )

    # Compute a mask whether the group exists.
    # (N, 9)
    group_exists = restype_rigidgroup_mask[aatype]

    # Compute a mask whether ground truth exists for the group
    # shape (N, 9, 3)
    gt_atoms_exist = np.take(all_atom_mask.astype(np.float32), base_atom_pos_idx,)

    gt_exists = np.min(gt_atoms_exist, axis=-1) * group_exists  # (N, 9)

    # Adapt backbone frame to old convention (mirror x-axis and z-axis).
    rots = np.tile(np.eye(3, dtype=np.float32), [_rc.num_frames, 1, 1])
    rots[0, 0, 0] = -1
    rots[0, 2, 2] = -1
    gt_frames = affine_mul_rots(gt_frames, rots)

    # The frames for ambiguous rigid groups are just rotated by 180 degree around
    # the x-axis. The ambiguous group is always the last chi-group.
    restype_rigidgroup_is_ambiguous = np.zeros(
        [_rc.full_num_residues, _rc.num_frames], dtype=np.float32
    )
    restype_rigidgroup_rots = np.tile(
        np.eye(3, dtype=np.float32), [_rc.full_num_residues, _rc.num_frames, 1, 1]
    )

    for resname, _ in _rc.residue_atom_renaming_swaps.items():
        restype = _rc.restype_order[_rc.restype_3to1[resname]]
        chi_idx = int(sum(_rc.chi_angles_mask[restype]) - 1)
        restype_rigidgroup_is_ambiguous[restype, chi_idx + 4] = 1
        restype_rigidgroup_rots[restype, chi_idx + 4, 1, 1] = -1
        restype_rigidgroup_rots[restype, chi_idx + 4, 2, 2] = -1

    # Gather the ambiguity information for each residue.
    residx_rigidgroup_is_ambiguous = restype_rigidgroup_is_ambiguous[aatype]
    residx_rigidgroup_ambiguity_rot = restype_rigidgroup_rots[aatype]

    # Create the alternative ground truth frames.
    alt_gt_frames = affine_mul_rots(gt_frames, residx_rigidgroup_ambiguity_rot)

    # reshape back to original residue layout
    gt_frames = np.reshape(gt_frames.numpy(), aatype_in_shape + (_rc.num_frames, 3, 4))
    gt_exists = np.reshape(gt_exists, aatype_in_shape + (_rc.num_frames,))
    group_exists = np.reshape(group_exists, aatype_in_shape + (_rc.num_frames,))
    residx_rigidgroup_is_ambiguous = np.reshape(
        residx_rigidgroup_is_ambiguous, aatype_in_shape + (_rc.num_frames,)
    )
    alt_gt_frames = np.reshape(
        alt_gt_frames.numpy(), aatype_in_shape + (_rc.num_frames, 3, 4,),
    )

    return {
        "rigidgroups_gt_frames": gt_frames,  # (..., 9, 3, 4)
        "rigidgroups_gt_exists": gt_exists,  # (..., 9)
        "rigidgroups_group_exists": group_exists,  # (..., 9)
        "rigidgroups_group_is_ambiguous": residx_rigidgroup_is_ambiguous,  # (..., 9)
        "rigidgroups_alt_gt_frames": alt_gt_frames,  # (..., 9, 3, 4)
    }


def atomf_to_torsion_angles(
    aatype: np.ndarray,  # (B, N)
    all_atom_positions: np.ndarray,  # (B, N, 65, 3)
    all_atom_mask: np.ndarray,  # (B, N, 65)
    placeholder_for_undefined=False,
) -> Dict[str, np.ndarray]:
    """Computes the 7 torsion angles (in sin, cos encoding) for each residue.
    The 7 torsion angles are in the order
    '[pre_omega, phi, psi, chi_1, chi_2, chi_3, chi_4]',
    here pre_omega denotes the omega torsion angle between the given amino acid
    and the previous amino acid.
    Args:
      aatype: Amino acid type, given as array with integers.
      all_atom_positions: atom37 representation of all atom coordinates.
      all_atom_mask: atom37 representation of mask on all atom coordinates.
      placeholder_for_undefined: flag denoting whether to set masked torsion
        angles to zero.
    Returns:
      Dict containing:
        * 'torsion_angles_sin_cos': Array with shape (B, N, 7, 2) where the final
          2 dimensions denote sin and cos respectively
        * 'alt_torsion_angles_sin_cos': same as 'torsion_angles_sin_cos', but
          with the angle shifted by pi for all chi angles affected by the naming
          ambiguities.
        * 'torsion_angles_mask': Mask for which chi angles are present.
    """

    # Map aatype > 28 to 'Unknown' (28).
    aatype = np.minimum(aatype, 28)
    prot_mask = aatype < _rc.num_prot

    # Compute the backbone angles.
    num_batch, num_res = aatype.shape

    pad = np.zeros([num_batch, 1, _rc.num_atoms, 3], np.float32)
    prev_all_atom_pos = np.concatenate([pad, all_atom_positions[:, :-1, :, :]], axis=1)

    pad = np.zeros([num_batch, 1, _rc.num_atoms], np.float32)
    prev_all_atom_mask = np.concatenate([pad, all_atom_mask[:, :-1, :]], axis=1)

    # For each torsion angle collect the 4 atom positions that define this angle.
    torsions_atom_pos = np.zeros(
        (num_batch, num_res, _rc.num_frames - 1, 4, 3), dtype=np.float32
    )
    torsion_angles_mask = np.zeros(
        (num_batch, num_res, _rc.num_frames - 1), dtype=np.float32
    )

    # shape (B, N, atoms=4, xyz=3)
    # Pre omega for proteins
    torsions_atom_pos[prot_mask, 0, :2] = prev_all_atom_pos[
        prot_mask, 1:3, :
    ]  # prev CA, C
    torsions_atom_pos[prot_mask, 0, 2:] = (
        all_atom_positions[prot_mask, :2, :],
    )  # this N, CA
    # Phi for proteins
    torsions_atom_pos[prot_mask, 1, :1] = prev_all_atom_pos[prot_mask, 2:3, :]  # prev C
    torsions_atom_pos[prot_mask, 1, 1:] = (
        all_atom_positions[prot_mask, :3, :],
    )  # this N, CA, C
    # Psi for proteins
    torsions_atom_pos[prot_mask, 2] = all_atom_positions[
        prot_mask, :4, :
    ]  # this N, CA, C, O

    # Collect the masks from these atoms.
    # Shape [batch, num_res]
    # Pre omega protein
    torsion_angles_mask[prot_mask, 0] = np.prod(
        prev_all_atom_mask[prot_mask, 1:3], axis=-1
    ) * np.prod(  # prev CA, C
        all_atom_mask[prot_mask, 0:2], axis=-1
    )  # this N, CA
    # Phi protein
    torsion_angles_mask[prot_mask, 1] = prev_all_atom_mask[
        prot_mask, 2
    ] * np.prod(  # prev C
        all_atom_mask[prot_mask, 0:3], axis=-1
    )  # this N, CA, C
    # Psi protein
    torsion_angles_mask[prot_mask, 2] = (
        np.prod(all_atom_mask[prot_mask, 0:3], axis=-1)
        * all_atom_mask[prot_mask, 4]  # this N, CA, C
    )  # this O

    # Grabbing atoms used for angles in NA is a bit complicated:
    na_aatype = aatype[~prot_mask]
    na_aatype -= _rc.num_prot  # Now 0 is DA, etc
    # Select atoms to compute angles. Shape: [batch, num_aa_res, angles=7, atoms=4].
    na_atom_indices = _rc.nuc_torsion_atom_indices[na_aatype] + np.arange(
        len(na_aatype) * _rc.num_atoms, step=_rc.num_atoms
    ).reshape(len(na_aatype), 1, 1)
    # Gather atom positions. Shape: [len(na_aatype), angles=7, atoms=4, xyz=3].
    torsions_atom_pos[~prot_mask, :7] = np.take(
        all_atom_positions[~prot_mask].reshape(-1, 3), na_atom_indices, axis=0
    )

    # Compute the mask in a similar way
    na_mask = _rc.nuc_torsion_atom_mask[na_aatype]
    na_atoms_mask = np.take(
        all_atom_mask[~prot_mask].reshape(-1), na_atom_indices, axis=0
    )
    na_atoms_mask = np.prod(na_atoms_mask, axis=-1)
    torsion_angles_mask[~prot_mask, :7] = na_mask * (na_atoms_mask).astype(np.float32)

    # Collect the atoms for the chi-angles.
    # Compute the table of chi angle indices. Shape: [restypes, chis=4, atoms=4].
    chi_atom_indices = _rc.chi_atom_indices
    # Select atoms to compute chis. Shape: [batch, num_res, chis=4, atoms=4].
    atom_indices = chi_atom_indices[aatype] + np.arange(
        num_res * num_batch * _rc.num_atoms, step=_rc.num_atoms
    ).reshape(num_batch, num_res, 1, 1)
    # Gather atom positions. Shape: [batch, num_res, chis=4, atoms=4, xyz=3].

    chis_atom_pos = np.take(all_atom_positions.reshape(-1, 3), atom_indices, axis=0)

    # Copy the chi angle mask, add the UNKNOWN residue. Shape: [restypes, 4].
    chi_angles_mask = list(_rc.chi_angles_mask)
    chi_angles_mask.append([0.0, 0.0, 0.0, 0.0])
    chi_angles_mask = np.asarray(chi_angles_mask)

    # Compute the chi angle mask. I.e. which chis angles exist according to the
    # aatype. Shape [batch, num_res, chis=4].
    chis_mask = chi_angles_mask[aatype]

    # Constrain the chis_mask to those chis, where the ground truth coordinates of
    # all defining four atoms are available.
    # Gather the chi angle atoms mask. Shape: [batch, num_res, chis=4, atoms=4].
    chi_angle_atoms_mask = np.take(all_atom_mask.reshape(-1), atom_indices, axis=0)
    # Check if all 4 chi angle atoms were set. Shape: [batch, num_res, chis=4].
    chi_angle_atoms_mask = np.prod(chi_angle_atoms_mask, axis=-1)
    chis_mask = chis_mask * (chi_angle_atoms_mask).astype(np.float32)

    torsions_atom_pos[prot_mask, -5:-1] = chis_atom_pos[prot_mask]
    torsions_atom_pos[~prot_mask, -1] = chis_atom_pos[~prot_mask][..., 0, :, :]

    torsion_angles_mask[prot_mask, -5:-1] = chis_mask[prot_mask]
    torsion_angles_mask[~prot_mask, -1] = chis_mask[~prot_mask][..., 0]

    # Create a frame from the first three atoms:
    # First atom: point on x-y-plane
    # Second atom: point on negative x-axis
    # Third atom: origin
    # Affine matrices (B, N, torsions=8, 3, 4)
    torsion_frames = affine_from_3_points(
        point_on_neg_x_axis=torch.Tensor(torsions_atom_pos[:, :, :, 1, :]),
        origin=torch.Tensor(torsions_atom_pos[:, :, :, 2, :]),
        point_on_xy_plane=torch.Tensor(torsions_atom_pos[:, :, :, 0, :]),
    )

    # Compute the position of the fourth atom in this frame (y and z coordinate
    # define the chi angle)
    # (B, N, torsions=7, 3)
    fourth_atom_rel_pos = affine_mul_vecs(
        invert_affine(torsion_frames), torch.Tensor(torsions_atom_pos[:, :, :, 3, :])
    ).numpy()

    # Normalize to have the sin and cos of the torsion angle.
    # np.ndarray (B, N, torsions=8, sincos=2)
    torsion_angles_sin_cos = np.stack(
        [fourth_atom_rel_pos[..., 2], fourth_atom_rel_pos[..., 1]], axis=-1
    )
    torsion_angles_sin_cos /= np.sqrt(
        np.sum(np.square(torsion_angles_sin_cos), axis=-1, keepdims=True) + 1e-8
    )

    # Mirror psi, because we computed it from the Oxygen-atom.
    torsion_angles_sin_cos[prot_mask] *= np.asarray(
        [1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    )[None, :, None]
    # torsion_angles_sin_cos[~prot_mask] *= np.asarray([1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])[
    #     None, :, None
    #  ]

    # Create alternative angles for ambiguous atom names.
    # Does not affect nucleotides
    chi_is_ambiguous = np.asarray(_rc.chi_pi_periodic)[aatype]

    mirror_torsion_angles = np.concatenate(
        [
            np.ones([num_batch, num_res, 3]),
            1.0 - 2.0 * chi_is_ambiguous,
            np.ones([num_batch, num_res, 1]),
        ],
        axis=-1,
    )
    alt_torsion_angles_sin_cos = (
        torsion_angles_sin_cos * mirror_torsion_angles[:, :, :, None]
    )

    if placeholder_for_undefined:
        # Add placeholder torsions in place of undefined torsion angles
        # (e.g. N-terminus pre-omega)
        placeholder_torsions = np.stack(
            [
                np.ones(torsion_angles_sin_cos.shape[:-1]),
                np.zeros(torsion_angles_sin_cos.shape[:-1]),
            ],
            axis=-1,
        )
        torsion_angles_sin_cos = torsion_angles_sin_cos * torsion_angles_mask[
            ..., None
        ] + placeholder_torsions * (1 - torsion_angles_mask[..., None])
        alt_torsion_angles_sin_cos = alt_torsion_angles_sin_cos * torsion_angles_mask[
            ..., None
        ] + placeholder_torsions * (1 - torsion_angles_mask[..., None])

    if num_batch == 1:
        torsion_angles_sin_cos = torsion_angles_sin_cos[0]
        alt_torsion_angles_sin_cos = alt_torsion_angles_sin_cos[0]
        torsion_angles_mask = torsion_angles_mask[0]

    return {
        "torsion_angles_sin_cos": torsion_angles_sin_cos,  # (B, N, 8, 2)
        "alt_torsion_angles_sin_cos": alt_torsion_angles_sin_cos,  # (B, N, 8, 2)
        "torsion_angles_mask": torsion_angles_mask,  # (B, N, 8)
    }


def torsion_angles_to_frames(
    aatype: np.ndarray,  # (N)
    backb_to_global: torch.Tensor,  # (N, 3, 4)
    torsion_angles_sin_cos: torch.Tensor,  # (N, 8, 2)
):  # (N, 9)
    """Compute rigid group frames from torsion angles.
    Jumper et al. (2021) Suppl. Alg. 24 "computeAllAtomCoordinates" lines 2-10
    Jumper et al. (2021) Suppl. Alg. 25 "makeRotX"
    Args:
      aatype: aatype for each residue
      backb_to_global: Rigid transformations describing transformation from
        backbone frame to global frame.
      torsion_angles_sin_cos: sin and cosine of the 8 torsion angles
    Returns:
      Frames corresponding to all the Sidechain Rigid Transforms
    """
    assert len(aatype.shape) == 1
    assert len(torsion_angles_sin_cos.shape) == 3
    assert torsion_angles_sin_cos.shape[1] == 8
    assert torsion_angles_sin_cos.shape[2] == 2

    device = torsion_angles_sin_cos.device
    # Gather the default frames for all rigid groups.
    # Affines with shape (N, 9, 3, 4)
    m = _rc.restype_rigid_group_default_frame[aatype]

    default_frames = affine_from_tensor4x4(torch.Tensor(m).to(device))

    # Create the rotation matrices according to the given angles (each frame is
    # defined such that its rotation is around the x-axis).
    sin_angles = torsion_angles_sin_cos[..., 0]
    cos_angles = torsion_angles_sin_cos[..., 1]

    # insert zero rotation for backbone group.
    (num_residues,) = aatype.shape
    sin_angles = torch.cat(
        [torch.zeros(num_residues, 1, device=device), sin_angles], dim=-1
    )
    cos_angles = torch.cat(
        [torch.ones(num_residues, 1, device=device), cos_angles], dim=-1
    )
    zeros = torch.zeros_like(sin_angles)
    ones = torch.ones_like(sin_angles)

    # all_rots are rotation_matrices with shape (N, 8, 3, 3)
    all_rots = fill_rotation_matrix(
        ones,
        zeros,
        zeros,
        zeros,
        cos_angles,
        -sin_angles,
        zeros,
        sin_angles,
        cos_angles,
    )

    # Apply rotations to the frames.
    all_frames = affine_mul_rots(default_frames, all_rots)

    # Nucleotide and protein torsions are defined very differently.
    # Thus, we need to split the computation
    all_frames_to_backb = torch.zeros_like(all_frames)
    prot_mask = aatype < _rc.num_prot

    # Proteins:
    # chi2, chi3, and chi4 frames do not transform to the backbone frame but to
    # the previous frame. So chain them up accordingly.

    all_frames_to_backb[prot_mask, :5] = all_frames[prot_mask, :5]
    all_frames_to_backb[prot_mask, 5] = affine_composition(
        all_frames_to_backb[prot_mask, 4], all_frames[prot_mask, 5]
    )  # chi2
    all_frames_to_backb[prot_mask, 6] = affine_composition(
        all_frames_to_backb[prot_mask, 5], all_frames[prot_mask, 6]
    )  # chi3
    all_frames_to_backb[prot_mask, 7] = affine_composition(
        all_frames_to_backb[prot_mask, 6], all_frames[prot_mask, 7]
    )  # chi4

    # Nucleotides:
    all_frames_to_backb[~prot_mask, :2] = all_frames[~prot_mask, :2]
    all_frames_to_backb[~prot_mask, 2] = affine_composition(
        all_frames_to_backb[~prot_mask, 1], all_frames[~prot_mask, 2]
    )  # beta
    all_frames_to_backb[~prot_mask, 3] = affine_composition(
        all_frames_to_backb[~prot_mask, 2], all_frames[~prot_mask, 3]
    )  # gamma
    all_frames_to_backb[~prot_mask, 4] = affine_composition(
        all_frames_to_backb[~prot_mask, 3], all_frames[~prot_mask, 4]
    )  # delta
    all_frames_to_backb[~prot_mask, 5] = affine_composition(
        all_frames_to_backb[~prot_mask, 3], all_frames[~prot_mask, 5]
    )  # nu2, from gamma
    all_frames_to_backb[~prot_mask, 6] = affine_composition(
        all_frames_to_backb[~prot_mask, 5], all_frames[~prot_mask, 6]
    )  # nu1
    all_frames_to_backb[~prot_mask, 7] = affine_composition(
        all_frames_to_backb[~prot_mask, 6], all_frames[~prot_mask, 7]
    )  # nu0
    all_frames_to_backb[~prot_mask, 8] = affine_composition(
        all_frames_to_backb[~prot_mask, 6], all_frames[~prot_mask, 8]
    )  # chi, from nu1

    # Create the global frames.
    # shape (N, 9, 3, 4)
    all_frames_to_global = affine_composition(
        backb_to_global[:, None], all_frames_to_backb
    )

    return all_frames_to_global


def frames_and_literature_positions_to_atomc_pos(
    aatype: np.ndarray, all_frames_to_global: torch.Tensor  # (N)  # (N, 9, 3, 4)
):  # (N, 23, 3)
    """Put atom literature positions (atom23 encoding) in each rigid group.
    Jumper et al. (2021) Suppl. Alg. 24 "computeAllAtomCoordinates" line 11
    Args:
      aatype: aatype for each residue.
      all_frames_to_global: All per residue coordinate frames.
    Returns:
      Positions of all atom coordinates in global frame.
    """

    device = all_frames_to_global.device
    # Pick the appropriate transform for every atom.
    residx_to_group_idx = _rc.restype_atomc_to_rigid_group[aatype]
    group_mask = torch.eye(9, device=device)[residx_to_group_idx.reshape(-1)].reshape(
        *residx_to_group_idx.shape, 9
    )  # shape (N, 23, 8)

    # Affines with shape (N, 23, 3, 4)
    # map_atoms_to_global = torch.sum(
    #     all_frames_to_global[:, None] * group_mask[..., None, None], dim=1
    # )
    map_atoms_to_global = torch.einsum(
        "nfij,npf->npij", all_frames_to_global, group_mask
    )

    # Gather the literature atom positions for each residue.
    # Vectors with shape (N, 23, 3)
    lit_positions = torch.Tensor(_rc.restype_atomc_rigid_group_positions[aatype]).to(
        device
    )

    # Transform each atom from its local frame to the global frame.
    # Vectors with shape (N, 23, 3)
    pred_positions = affine_mul_vecs(map_atoms_to_global, lit_positions)

    # Mask out non-existing atoms.
    mask = torch.Tensor(_rc.restype_atomc_mask[aatype]).to(device)
    pred_positions = pred_positions * mask[..., None]

    return pred_positions


def frames_and_literature_positions_to_atom3_pos(
    aatype: np.ndarray, all_frames_to_global: torch.Tensor  # (N)  # (N, 3, 4)
):  # (N, 3, 3)
    """Put atom literature positions (atom3 encoding) in each rigid group.
    It should be in N,CA,C order.
    Similar to Jumper et al. (2021) Suppl. Alg. 24 "computeAllAtomCoordinates" line 11
    Args:
      aatype: aatype for each residue.
      all_frames_to_global: All per residue coordinate frames.
    Returns:
      Positions of all atom coordinates in global frame.
    """
    if isinstance(aatype, torch.Tensor):
        aatype = aatype.cpu().detach().numpy()

    device = all_frames_to_global.device

    # Gather the literature atom positions for each residue.
    # Vectors with shape (N, 3, 3)
    lit_positions = torch.Tensor(_rc.restype_atom3_rigid_group_positions[aatype]).to(
        device
    )

    # Transform each atom from its local frame to the global frame.
    # Vectors with shape (N, 3, 3)
    pred_positions = affine_mul_vecs(all_frames_to_global, lit_positions)
    return pred_positions


def add_lm_embeddings_to_protein(
    input_protein: Protein, lm_embeddings: np.ndarray
) -> Protein:
    assert len(lm_embeddings) == input_protein.unified_seq_len
    protein_dict_without_lm = dict(
        [
            (k, v)
            for (k, v) in input_protein.__dict__.items()
            if "residue_to_lm_embedding" not in k
        ]
    )
    new_protein = Protein(
        **protein_dict_without_lm, residue_to_lm_embedding=lm_embeddings,
    )
    return new_protein


def load_protein_from_prot(file_path):
    with open(file_path, "rb") as f:
        prot = pickle.load(f)
    return Protein(**prot)


def dump_protein_to_prot(input_protein, file_path):
    with open(file_path, "wb") as f:
        pickle.dump(input_protein.__dict__, f)


def get_sequence_context_from_idx(idx_arr, num_residues, residue_to_seq_id, context=20):
    sorted_idx = np.argsort(idx_arr)
    idx_to_seq_context = np.zeros_like(idx_arr, dtype=np.int64)
    sequence_context = list(
        range(max(idx_arr[sorted_idx[0]] - context, 0), idx_arr[sorted_idx[0]])
    )
    last_node = idx_arr[sorted_idx[0]]
    idx_to_seq_context[sorted_idx[0]] = len(sequence_context)
    for i, idx in enumerate(idx_arr[sorted_idx[1:]], start=1):
        if idx - last_node < context:
            sequence_context.extend(range(last_node, idx))
        else:
            sequence_context.extend(
                range(last_node, min(last_node + context, num_residues))
            )
            sequence_context.extend(
                range(max(idx - context, sequence_context[-1]), idx)
            )
        idx_to_seq_context[sorted_idx[i]] = len(sequence_context)
        last_node = idx
    sequence_context.extend(range(last_node, min(last_node + context, num_residues)))
    non_unique_seq = residue_to_seq_id[np.array(sequence_context, dtype=np.int64)]
    unique_seq, reverse_idx = np.unique(non_unique_seq, return_inverse=True)
    return unique_seq, reverse_idx[idx_to_seq_context]


def make_up_fr_from_prot(prot: Protein):
    final_results = {}
    final_results["existence_mask"] = np.ones(len(prot.aatype)) > 0.5
    final_results["aa_logits"] = np.zeros((len(prot.aatype, 28)))
    final_results["aa_logits"][np.arange(len(prot.aatype)), prot.aatype] = 2
    final_results["pred_affines"] = prot.rigidgroups_gt_frames[:, 0]


def slice_protein(protein: Protein, slice_array: np.ndarray) -> Protein:
    num_res = len(protein.aatype)
    for key in PROTEIN_KEYS:
        value = getattr(protein, key)
        if hasattr(value, "shape") and value.shape[0] == num_res:
            setattr(protein, key, value[slice_array])
    return protein


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Generate .prot files from PDB/mmCIF")
    parser.add_argument("--input", required=True, help="The path to the PDB/mmCIF file")
    parser.add_argument(
        "--output",
        required=True,
        help="The path to save the .prot file (could end with anything)",
    )
    args = parser.parse_args()

    protein = get_protein_from_file_path(args.input)
    dump_protein_to_prot(protein, args.output)

    print(f"Job successful, output at {args.output}")
