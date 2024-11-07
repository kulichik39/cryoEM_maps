import os
import numpy as np
import mrcfile
from subprocess import PIPE, Popen


def run_script_inside_chimera(
    script_name,
    script_path=os.getcwd() + os.path.sep + "chimera_scripts",
    option_names=(),
    option_values=(),
):
    """
    Runs given python script inside Chimera through the subprocess (simulates running
    script from the terminal command:
    chimera --nogui --nostatus --script "name_of_script.py").
    Returns an object of the corresponding subprocess.

    Params:
    script_name - name of the python script
    script_path - path to the python script (excluding its name)
    option_names - names of the options provided to the script
    option_values - values for the options provided to the script

    Returns:
    p - object of the Popen class corresponding to the script's subprocess
    """

    assert script_name.endswith(
        ".py"
    ), "The name of the python script must end with .py"

    assert len(option_names) == len(
        option_values
    ), "Mismatch between number of options and number of values provided"

    # full path to the python script (including its name)
    path_full = script_path + os.path.sep + script_name

    # construct the line of the terminal command that describes the script
    script_line = path_full
    for i in range(len(option_names)):
        script_line += f" {option_names[i]}" + f" {option_values[i]}"

    p = Popen(
        [
            "chimera",
            "--nogui",
            "--nostatus",
            "--script",
            script_line,
        ],
        stdout=PIPE,
        stderr=PIPE,
    )

    return p


def prepare_raw_pdb_file(
    raw_pdb_filename,
    raw_molecule_path=os.getcwd() + os.path.sep + "raw_molecule_data",
    molecule_path=os.getcwd() + os.path.sep + "molecule_data",
):
    """
    Prepares a raw .pdb molecule file to feed into Chimera: deletes all unnecessary rows.

    Params:
    raw_pdb_filename - name of the raw input .pdb file
    raw_molecule_path - path to the directory where raw molecule files are stored
    molecule_path - path to the directory where processed (prepared) molecule files are stored

    Returns:
    prepared_pdb_filename - name of the prepared_pdb_file
    """

    assert raw_pdb_filename.endswith(".pdb"), "pdb filename must end with .pdb!"

    # create directory for prepared molecule files if it doesn't exist
    if not os.path.exists(molecule_path):
        os.mkdir(molecule_path)

    # construct full path to the input raw pdb file
    raw_pdb_path_full = raw_molecule_path + os.path.sep + raw_pdb_filename

    # construct full path to the prepared pdb file
    prepared_pdb_filename = "prepared_" + raw_pdb_filename
    pdb_path_full = molecule_path + os.path.sep + prepared_pdb_filename

    with open(raw_pdb_path_full, "r") as raw_file:
        with open(pdb_path_full, "w") as prepared_file:
            for raw_line in raw_file:
                # Keep only the lines corresponding to atoms
                if raw_line.startswith("ATOM"):
                    prepared_file.write(raw_line)

    return prepared_pdb_filename


def read_density_data_mrc(
    density_filename,
    density_path=os.getcwd() + os.path.sep + "density_maps",
):
    """
    Reads density map data from the given .mrc file and converts it to numpy array.

    Params:
    density_filename - name of the file with density data
    density_path - path to the directory where density file is stored

    Returns:
    density - numpy array with density data
    """

    assert density_filename.endswith(".mrc"), "mrc filename must end with .mrc!"

    # full path to the density file (including its name)
    path_full = density_path + os.path.sep + density_filename

    with mrcfile.open(path_full) as mrc:
        density = mrc.data

    return density


def read_density_data_txt(
    density_filename,
    density_path=os.getcwd() + os.path.sep + "density_maps",
):
    """
    Reads density map data from the given .txt file and converts it to numpy array.

    Params:
    density_filename - name of the file with density data
    density_path - path to the directory where density file is stored

    Returns:
    density - numpy array with density data
    """

    assert density_filename.endswith(".txt"), "txt filename must end with .txt!"

    # full path to the density file (including its name)
    path_full = density_path + os.path.sep + density_filename

    # fetch dimensions of the density matrix from the file
    dims = tuple()
    with open(path_full, "r") as file:
        first_line = file.readline()
        first_line = first_line.split(":")[-1]  # dimensions are written here
        first_line = first_line[1:-2]  # remove \n and parentheses
        dims = tuple(map(int, first_line.split(", ")))

    density = np.loadtxt(path_full)
    density = density.reshape(dims)

    return density


def delete_extension_from_filename(filename):
    """
    Deletes extension (symbols after .) from the given filename.

    Params:
    filename - name of the given file

    Returns:
    name of the file without extension
    """

    return ".".join(filename.split(".")[:-1])
