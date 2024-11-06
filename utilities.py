import os
import numpy as np
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


def read_density_data(
    density_filename,
    density_path=os.getcwd() + os.path.sep + "density_maps",
):
    """
    Reads density map data from the given file and converts it to numpy array.

    Params:
    density_filename - name of the file with density data
    density_path - path to the directory where density file is stored

    Returns:
    density - numpy array with density data
    """

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
