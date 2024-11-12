import os
import numpy as np
from utilities import (
    read_molecule,
    delete_extension_from_filename,
    save_all_conformers_to_pdb,
    extract_filename_from_full_path,
    read_density_data_mrc,
)

from internal_flexibility import (
    generate_conformers,
    compute_density_for_multiple_conformers,
)

from missing_parts import random_delete_atoms_from_multiple_pdb_files

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

    # randomly delete atoms from the generated conformers
    delete_prob = 0.2  # probability of atom deleting
    # path to store cofomormers with deleted atoms
    delatoms_conformer_path = os.getcwd() + os.path.sep + "delatoms_conformers_data"
    delatoms_conformer_path_list = random_delete_atoms_from_multiple_pdb_files(
        conformer_path_list,
        delatoms_molecule_path=delatoms_conformer_path,
        delete_prob=delete_prob,
    )

    # compute density maps for generated conformers
    print(delatoms_conformer_path_list)
    p_list = compute_density_for_multiple_conformers(delatoms_conformer_path_list)

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
                    extract_filename_from_full_path(delatoms_conformer_path_list[i])
                )
                + ".mrc"
            )
        density = read_density_data_mrc(density_filename)
        print(density.shape)
        print(np.nonzero(density))
