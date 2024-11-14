import os
import numpy as np
from utilities import (
    compute_density_map_in_chimera,
    read_density_data_mrc,
    delete_extension_from_filename,
    extract_filename_from_full_path,
)


def random_delete_atoms_from_pdb_file(
    pdb_path_full,
    delatoms_molecule_path=os.getcwd() + os.path.sep + "delatoms_molecule_data",
    delete_prob=0.2,
):
    """
    Randomly deletes atoms from the given pdb file to simulate the case when
    some parts of the molecule are completely missed.

    Params:
    pdb_path_full - full path to the input .pdb file (icnluding its name)
    delatoms_molecule_path - path to the directory where files with some deleted atoms are
    stored
    delete_prob - probability for deleting an atom

    Returns:
    delatoms_pdb_path_full - full path to the new .pdb file with some atoms deleted
    """

    # exctract pdf filename from the input full path
    pdb_filename = extract_filename_from_full_path(pdb_path_full)

    assert pdb_filename.endswith(".pdb"), "pdb filename must end with .pdb!"

    # create directory for the molecule files with some deleted atoms if it doesn't exist
    if not os.path.exists(delatoms_molecule_path):
        os.mkdir(delatoms_molecule_path)

    # construct full path to the pdb file with some deleted atoms
    delatoms_pdb_path_full = delatoms_molecule_path + os.path.sep + pdb_filename

    with open(pdb_path_full, "r") as input_file:
        with open(delatoms_pdb_path_full, "w") as output_file:
            for input_line in input_file:

                if input_line.startswith("ATOM") or input_line.startswith("HETATM"):
                    # if the current line corresponds to an atom, decide whether this atom should
                    # be deleted by using uniform random number between 0 and 1.
                    random_val = np.random.rand()
                    if random_val > delete_prob:
                        output_file.write(input_line)
                else:
                    output_file.write(input_line)

    return delatoms_pdb_path_full


def random_delete_atoms_from_multiple_pdb_files(
    pdb_path_full_list,
    delatoms_molecule_path=os.getcwd() + os.path.sep + "delatoms_molecule_data",
    delete_prob=0.2,
):
    """
    Randomly deletes atoms from the several given pdb files to simulate the case when
    some parts of the molecules are completely missed.

    Params:
    pdb_path_full_list - list with full paths to the input .pdb files (icnluding their names)
    delatoms_molecule_path - path to the directory where files with some deleted atoms are
    stored
    delete_prob - probability for deleting an atom

    Returns:
    delatoms_pdb_path_full_list - list with full paths to the new .pdb files with some atoms
    deleted
    """

    delatoms_pdb_path_full_list = []  # list to store full paths of the output pdb files

    for pdb_path_full in pdb_path_full_list:
        delatoms_pdb_path_full = random_delete_atoms_from_pdb_file(
            pdb_path_full,
            delatoms_molecule_path=delatoms_molecule_path,
            delete_prob=delete_prob,
        )
        delatoms_pdb_path_full_list.append(delatoms_pdb_path_full)

    return delatoms_pdb_path_full_list


if __name__ == "__main__":

    # construct input paths
    input_filename = "1c3b_ligand.pdb"
    input_path = os.getcwd() + os.path.sep + "raw_molecule_data"
    input_path_full = input_path + os.path.sep + input_filename

    # delete some atoms from the input pdb file and write a new file
    del_pdb_path = random_delete_atoms_from_pdb_file(
        input_path_full
    )  # returns full path to new pdb file

    # run python script inside Chimera to compute density map
    density_resolution = 3.5  # resolution of the density map (in Angstrom)
    p = compute_density_map_in_chimera(
        del_pdb_path, density_resolution=density_resolution
    )
    stdout, stderr = p.communicate()  # output from the subprocess

    # if the subprocess finished with errors, check logs in chimera_logs folder
    if stderr:
        print(stderr)
        print("Failed to compute density map. Check logs in chimera_logs folder.")

    # if the subprocess finished correctly, check obtained denisty map
    # in the densisty_maps folder
    else:
        density_filename = (
            delete_extension_from_filename(
                extract_filename_from_full_path(del_pdb_path)
            )
            + ".mrc"
        )
        density = read_density_data_mrc(density_filename)
        print(density.shape)
        print(np.nonzero(density))
