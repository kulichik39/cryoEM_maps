import os
import numpy as np

from utilities import (
    run_script_inside_chimera,
    read_density_data_mrc,
    prepare_raw_pdb_file,
    delete_extension_from_filename,
)


def random_delete_atoms_from_pdb_file(
    pdb_filename,
    molecule_path=os.getcwd() + os.path.sep + "molecule_data",
    delete_prob=0.2,
):
    """
    Randomly deletes atoms from the given pdb file to simulate the case when
    some part of the molecule are completely missed.

    Params:
    pdb_filename - name of the input .pdb file
    molecule_path - path to the directory where processed (prepared) molecule files are stored.
    check prepare_raw... functions in utilities.py for more info
    delete_prob - probability for deleting an atom

    Returns:
    pdb_path_full_modified - full path to the new .pdb file with some atoms deleted
    """

    assert pdb_filename.endswith(".pdb"), "pdb filename must end with .pdb!"

    # construct full path to the input pdb file
    pdb_path_full = molecule_path + os.path.sep + pdb_filename

    # construct full path to the modified pdb file
    pdb_path_full_modified = molecule_path + os.path.sep + "del_atoms_" + pdb_filename

    with open(pdb_path_full, "r") as input_file:
        with open(pdb_path_full_modified, "w") as output_file:
            for input_line in input_file:
                # decide whether an atom should be deleted by using uniform random
                # number between 0 and 1.
                random_val = np.random.rand()
                if random_val > delete_prob:
                    output_file.write(input_line)

    return pdb_path_full_modified


if __name__ == "__main__":

    raw_input_pdb = "1c3b_ligand.pdb"  # name of the input pdb file

    # prepare raw pdb file to feed into Chimera, check prepare_raw_pdb_file() function
    # in utilities.py for more details
    prepared_input_pdb = prepare_raw_pdb_file(raw_input_pdb)

    # delete some parts of the prepared input pdb file and write a new file
    del_pdb_path = random_delete_atoms_from_pdb_file(
        prepared_input_pdb
    )  # full path to new pdb file

    # run python script inside Chimera
    scipt_name = (
        "chimera_density_map.py"  # name of the python script to run inside Chimera
    )
    # options for the python script
    option_names = ("-i",)
    option_values = (del_pdb_path,)

    # run python script inside Chimera as a subprocess
    # store density maps in the density_maps folder
    p = run_script_inside_chimera(
        scipt_name, option_names=option_names, option_values=option_values
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
            "density_map"
            + "_"
            + delete_extension_from_filename(del_pdb_path.split(os.path.sep)[-1])
            + ".mrc"
        )
        density = read_density_data_mrc(density_filename)
        print(density.shape)
        print(np.nonzero(density))
