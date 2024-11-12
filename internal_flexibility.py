import os
import numpy as np
from utilities import (
    read_molecule,
    save_all_conformers_to_pdb,
    read_density_data_mrc,
    compute_density_map_in_chimera,
    delete_extension_from_filename,
    extract_filename_from_full_path,
)
from rdkit import Chem
from rdkit.Chem import AllChem


def generate_conformers(
    mol, n_confs, add_Hs=True, use_small_ring_torsions=True, prune_rms_tresh=1.0
):
    """
    Generates n_confs conformers for the give molecule using RDKit library.

    Params:
    mol - input RDKit molecule object
    n_confs - number of conformers to generate
    add_Hs - whether to add hydrogen atoms to the given molecule before conformers generation
    use_small_ring_torsions - whether to include additional small ring torsion potentials
    prune_rms_tresh - threshold value for RMSD pruning of the conformers

    Returns:
    mol - RDKit molecule object with generated conformers
    """

    if add_Hs:  # add hydrogen atoms if needed
        mol = Chem.AddHs(mol)

    # specify parameters for the conformers generation
    params = AllChem.ETKDGv3()
    params.useSmallRingTorsions = use_small_ring_torsions
    params.pruneRMsThresh = prune_rms_tresh

    AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)

    return mol


def compute_density_for_one_conformer(
    conformer_path_full,
    script_path=os.getcwd() + os.path.sep + "chimera_scripts",
):
    """
    Computes density map for the one specified conformer by running python
    script inside Chimera.

    Params:
    conformer_path_full - full path to the input conformer file (including its name)
    script_path - path to the folder where the python script for density computation is stored
    (excluding its name)

    Returns:
    p - object of the Popen class corresponding to the Chimera's subprocess
    """

    # run scipt inside Chimera and return the subprocess object
    p = compute_density_map_in_chimera(conformer_path_full, script_path=script_path)

    return p


def compute_density_for_multiple_conformers(
    conformer_path_full_list,
    script_path=os.getcwd() + os.path.sep + "chimera_scripts",
):
    """
    Computes density maps for several specified conformers by running python
    script inside Chimera.

    Params:
    conformer_path_full_list - list with full paths to the input conformer files
    (including their names)
    script_path - path to the folder where the python script for density computation is stored
    (excluding its name)

    Returns:
    p_list - list of Popen objects corresponding to Chimera's subprocesses
    """

    p_list = []

    # iterate over the given input conformers files
    for conformer_path_full in conformer_path_full_list:
        p = compute_density_for_one_conformer(
            conformer_path_full, script_path=script_path
        )
        p_list.append(p)

    return p_list


def compute_average_map():
    pass


if __name__ == "__main__":

    # construct input paths
    input_filename = "5SD5_HWI_ligand.sdf"
    input_path = os.getcwd() + os.path.sep + "raw_molecule_data"

    # read molecule data using RDKit
    mol = read_molecule(input_filename, file_path=input_path, remove_Hs=True)

    # generate conformers
    n_confs = 5
    mol = generate_conformers(
        mol,
        n_confs=n_confs,
        add_Hs=True,
        use_small_ring_torsions=True,
        prune_rms_tresh=1.0,
    )

    # save generated conformers to files
    conformers_path = os.getcwd() + os.path.sep + "raw_conformers_data"
    base_conformer_filename = delete_extension_from_filename(input_filename) + ".pdb"
    conformer_path_list = save_all_conformers_to_pdb(
        mol, base_conformer_filename, pdb_path=conformers_path
    )

    # compute density maps for generated conformers
    p_list = compute_density_for_multiple_conformers(conformer_path_list)

    # check computed density maps
    for i in range(len(p_list)):
        stdout, stderr = p_list[i].communicate()

        # if the subprocess finished with errors, check logs in chimera_logs folder
        if stderr:
            print(stderr)
            print("Failed to compute density map. Check logs in chimera_logs folder.")

        # if the subprocess finished correctly, check obtained denisty map
        # in the densisty_maps folder
        else:
            density_filename = (
                delete_extension_from_filename(
                    extract_filename_from_full_path(conformer_path_list[i])
                )
                + ".mrc"
            )
        density = read_density_data_mrc(density_filename)
        print(density.shape)
        print(np.nonzero(density))
