import os
import numpy as np
from utilities import (
    read_molecule,
    delete_extension_from_filename,
    save_all_conformers_to_pdb,
    extract_filename_from_full_path,
    read_density_data_mrc,
    group_conformers_to_single_file,
    compute_density_map_in_chimera,
)

from internal_flexibility import (
    generate_conformers,
)

from missing_parts import random_delete_atoms_from_pdb_file

if __name__ == "__main__":

    # construct input paths
    input_filename = "1c3b_ligand.pdb"
    input_path = os.getcwd() + os.path.sep + "raw_molecule_data"

    # read molecule data using RDKit
    mol = read_molecule(input_filename, file_path=input_path, remove_Hs=True)

    # generate conformers
    n_confs = 5
    mol = generate_conformers(
        mol,
        n_confs=n_confs,
        add_Hs=False,
        use_small_ring_torsions=True,
        prune_rms_tresh=1.0,
    )

    # save generated conformers to files
    conformers_path = os.getcwd() + os.path.sep + "raw_conformers_data"
    base_conformer_filename = delete_extension_from_filename(input_filename) + ".pdb"
    conformer_path_list = save_all_conformers_to_pdb(
        mol, base_conformer_filename, pdb_path=conformers_path
    )

    # group conformers' files into a single file to generate an "average" map
    group_filename = (
        "group_conformers_" + delete_extension_from_filename(input_filename) + ".pdb"
    )
    group_conformers_path = group_conformers_to_single_file(
        conformer_path_list,
        group_filename,
        group_path=conformers_path,
        delete_input=True,
    )

    # randomly delete atoms from the generated conformers
    delete_prob = 0.2  # probability of atom deleting
    # path to the folder where conformers with deleted atoms are stored
    delatoms_conformer_folder = os.getcwd() + os.path.sep + "delatoms_conformers_data"
    delatoms_conformer_path = random_delete_atoms_from_pdb_file(
        group_conformers_path,
        delatoms_molecule_path=delatoms_conformer_folder,
        delete_prob=delete_prob,
    )

    # compute "average" density map for the grouped conformers
    density_resolution = 3.5  # resolution of the density map (in Angstrom)
    p = compute_density_map_in_chimera(
        delatoms_conformer_path, density_resolution=density_resolution
    )

    # check the density
    stdout, stderr = p.communicate()

    # if the subprocess finished with errors, check logs in chimera_logs folder
    if stderr:
        print(stderr)
        print("Failed to compute density map. Check logs in chimera_logs folder.")

    # if the subprocess finished correctly, check obtained denisty map
    # in the densisty_maps folder
    else:
        density_filename = (
            delete_extension_from_filename(
                extract_filename_from_full_path(group_conformers_path)
            )
            + ".mrc"
        )
        density = read_density_data_mrc(density_filename)
        print(density.shape)
        print(np.nonzero(density))
