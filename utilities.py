import os
import numpy as np
import mrcfile
from subprocess import Popen
from rdkit import Chem
from datetime import datetime, timezone


def run_script_inside_chimera(
    script_name,
    script_path=os.getcwd() + os.path.sep + "chimera_scripts",
    option_names=(),
    option_values=(),
    stderr_file=None
):
    """
    Runs given python script inside Chimera through a subprocess (simulates running
    script from the terminal command:
    chimera --nogui --nostatus --script "name_of_script.py").
    NOTE: if some option in option_names doesn't have a value, an empty string ""
    must be provided to the option_values.
    Returns an object of the corresponding subprocess.

    Params:
    script_name - name of the python script
    script_path - path to the python script (excluding its name)
    option_names - names of the options provided to the script
    option_values - values for the options provided to the script
    stderr_file - file object to log errors from the subprocess if any

    Returns:
    p - object of the Popen class corresponding to the Chimera's subprocess
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
        if option_values[i]:  # if option value is not an empty string
            option = f" {option_names[i]}" + f" {option_values[i]}"
        else:
            option = f" {option_names[i]}"
        script_line += option

    p = Popen(
        [
            "chimera",
            "--nogui",
            "--nostatus",
            "--script",
            script_line,
        ],
        stderr=stderr_file,
    )

    return p


def compute_density_map_in_chimera(
    molecule_path_full,
    density_path_full,
    is_log=False,
    log_path=os.getcwd() + os.path.sep + "chimera_logs",
    density_resolution=1.0,
    script_path=os.getcwd() + os.path.sep + "chimera_scripts",
    stderr_file=None
):
    """
    Computes density map for the given molecule/conformer using Chimera.
    To achieve this, runs python script "chimera_density_map.py" inside Chimera through a
    subprocess.

    Params:
    molecule_path_full - full path to the input molecule file (including its name)
    density_path_full - full path to the output mrc density file (including its name)
    is_log - should we write logs for Chimera scripts
    log_path - path to the folder where the log file will be stored (excluding the file's name which will be created automatically)
    density_resolution - desired resolution of the map (in Angstrom)
    script_path - path to the folder with python script (excluding its name)
    stderr_file - file object to log errors from the subprocess if any

    Returns:
    p - object of the Popen class corresponding to the Chimera's subprocess
    """

    # name of the python script to run inside Chimera
    scipt_name = "chimera_density_map.py"

    # options for the python script
    option_names = ["-i", "-r", "-o"]
    option_values = [molecule_path_full, density_resolution, density_path_full]
    if is_log:
        option_names.append("-l")
        option_values.append(log_path)

    # returns object corresponding to the subprocess
    p = run_script_inside_chimera(
        scipt_name,
        script_path=script_path,
        option_names=option_names,
        option_values=option_values,
        stderr_file=stderr_file
    )

    return p


# NOTE: for now, the function commented out below shouldn't be used in the code

# def prepare_raw_pdb_file(
#     raw_pdb_filename,
#     raw_molecule_path=os.getcwd() + os.path.sep + "raw_molecule_data",
#     molecule_path=os.getcwd() + os.path.sep + "molecule_data",
# ):
#     """
#     Prepares a raw .pdb molecule file to feed into Chimera: deletes all unnecessary rows.

#     Params:
#     raw_pdb_filename - name of the raw input .pdb file
#     raw_molecule_path - path to the directory where raw molecule files are stored
#     molecule_path - path to the directory where processed (prepared) molecule files are stored

#     Returns:
#     prepared_pdb_filename - name of the prepared_pdb_file
#     """

#     assert raw_pdb_filename.endswith(".pdb"), "pdb filename must end with .pdb!"

#     # create directory for prepared molecule files if it doesn't exist
#     create_folder(molecule_path)

#     # construct full path to the input raw pdb file
#     raw_pdb_path_full = raw_molecule_path + os.path.sep + raw_pdb_filename

#     # construct full path to the prepared pdb file
#     prepared_pdb_filename = "prepared_" + raw_pdb_filename
#     pdb_path_full = molecule_path + os.path.sep + prepared_pdb_filename

#     with open(raw_pdb_path_full, "r") as raw_file:
#         with open(pdb_path_full, "w") as prepared_file:
#             for raw_line in raw_file:
#                 # Keep only the lines corresponding to atoms
#                 if raw_line.startswith("ATOM") or raw_line.startswith("HETATM"):
#                     prepared_file.write(raw_line)

#     return prepared_pdb_filename


def read_density_data_mrc(
    density_path_full
):
    """
    Reads density map data from the given .mrc file and converts it to numpy array.

    Params:
    denisty_path_full - full path to the density file (including its name)

    Returns:
    density - numpy array with density data
    """

    assert density_path_full.endswith(".mrc"), "mrc filename must end with .mrc!"

    with mrcfile.open(density_path_full) as mrc:
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


def extract_filename_from_full_path(path_full):
    """
    Extracts filename from the full (including the name) path to file.

    Params:
    full_path - full path to the file

    Returns:
    filename extracted from the full path
    """

    return path_full.split(os.path.sep)[-1]


def extract_format_from_filename(filename):
    """
    Gets file's format from the given filename.

    Params:
    filename - name of the file

    Returns:
    format of the file
    """

    return filename.split(".")[-1]


def read_molecules_from_sdf(
    sdf_path_full,
    n_mols=1,
    remove_Hs=True,
):
    """
    Reads molecules data from the given .sdf file using RDKit library.

    Params:
    sdf_path_full - full path to the input .sdf file (including its name)
    n_mols - expected number of molecules to read from the file
    remove_Hs - whether to remove hydrogen atoms from the molecule when reading

    Returns:
    mols - list with RDKit molecule objects obtained from the given file
    """

    assert sdf_path_full.endswith(".sdf"), "sdf filename must end with .sdf!"

    mols = []  # list to store RDKit molecule objects

    with Chem.SDMolSupplier(sdf_path_full, removeHs=remove_Hs) as suppl:
        for mol in suppl:
            if mol is None:
                print(
                    f"Failed to read one of the molecules from the sdf file in: {sdf_path_full}"
                )
                continue
            mols.append(mol)

    assert (
        len(mols) == n_mols
    ), f"Expected {n_mols} molecules from the sdf file, but read {len(mols)} molecules. File path: {sdf_path_full}"

    return mols


def read_molecule_from_pdb(
    pdb_path_full,
    remove_Hs=True,
):
    """
    Reads one molecule data from the given .pdb file using RDKit library.
    NOTE: even if you have several molecules in the given file, RDKit will
    read only the first one.

    Params:
    pdb_path_full - full path to the input .pdb file (including its name)
    remove_Hs - whether to remove hydrogen atoms from the molecule when reading

    Returns:
    mol - RDKit molecule object obtained from the given file
    """

    assert pdb_path_full.endswith(".pdb"), "pdb filename must end with .pdb!"

    mol = Chem.MolFromPDBFile(pdb_path_full, removeHs=remove_Hs)

    return mol


def read_molecule(
    path_full, remove_Hs=True
):
    """
    Reads one molecule data from the given file and converts it to an RDKit molecule object.

    Params:
    path_full - full path to the input molecule file
    remove_Hs - whether to remove hydrogen atoms from the molecule when reading

    Returns:
    mol - RDKit molecule object obtained from the given file
    """

    # exctract name of the file from full path
    filename = extract_filename_from_full_path(path_full)

    # extract format of the name
    file_format = extract_format_from_filename(filename)

    match file_format:
        case "sdf":
            mols = read_molecules_from_sdf(
                path_full, n_mols=1, remove_Hs=remove_Hs
            )
            return mols[0]

        case "pdb":
            mol = read_molecule_from_pdb(
                path_full, remove_Hs=remove_Hs
            )
            return mol

        case default:
            raise RuntimeError(
                f"Reading from .{file_format} molecule files is not implemented!"
            )


def save_one_conformer_to_pdb(
    mol,
    conf_id,
    pdb_filename,
    pdb_path=os.getcwd() + os.path.sep + "raw_conformers_data",
):
    """
    Saves specified (by conf_id) conformer of the given RDKit molecule to a .pdb file.

    Params:
    mol - RDKit molecule object
    conf_id - id of the molecule's conformer to save. -1 indicates the last conformer
    pdb_filename - name of the pdb file
    pdb_path - path to the folder where the pdb file should stored (excluding its name)

    Returns:
    pdb_path_full - full path to the conformer's pdb file (including its name)
    """

    assert pdb_filename.endswith(".pdb"), "pdb filename must end with .pdb!"

    # create the folder for saving if it doesn't exist
    create_folder(pdb_path)

    # construct full path to the pdb file (including its name)
    pdb_path_full = pdb_path + os.path.sep + pdb_filename

    Chem.MolToPDBFile(
        mol,
        pdb_path_full,
        confId=conf_id,
    )

    return pdb_path_full


def save_all_conformers_to_pdb(
    mol,
    base_pdb_filename,
    pdb_path=os.getcwd() + os.path.sep + "raw_conformers_data",
):
    """
    Saves all conformers of the given RDKit molecule to .pdb files.
    Each conformer is saved in a separate file.

    Params:
    mol - RDKit molecule object
    base_pdb_filename - base name of the output pdb files. "conformer_{conformer_id}" will
    be added to this name
    pdb_path - path to the folder where the pdb files should be stored

    Returns:
    pdb_path_full_list - list with full paths to the conformers' pdb files
    """

    assert base_pdb_filename.endswith(".pdb"), "pdb filename must end with .pdb!"

    pdb_path_full_list = []  # list with full paths to the output conformers' files

    # iterate over all available conformers
    conf_ids = [x.GetId() for x in mol.GetConformers()]
    for conf_id in conf_ids:
        # construct the name of the pdb file for the given conformer
        pdb_filename = f"conf_{conf_id}_" + base_pdb_filename
        pdb_path_full = save_one_conformer_to_pdb(
            mol, conf_id, pdb_filename, pdb_path=pdb_path
        )
        pdb_path_full_list.append(pdb_path_full)

    return pdb_path_full_list


def group_conformers_to_single_file(
    path_full_list,
    group_filename,
    group_path=os.getcwd() + os.path.sep + "raw_conformers_data",
    delete_input=True,
):
    """
    Groups multiple conformers data into a single file. The paths to the input conformers files are
    provided in the path_full_list.

    Params:
    path_full_list - list with full paths to the input conformers' files to group (including their names)
    group_filename - name of the file where the conformers will be written
    group_path - path to the group file (excluding its name)
    delete_input - whether to delete input conformers' files

    Returns:
    group_path_full - full path to the file where conformers are grouped (including its name)
    """

    # exctract format of the group file
    group_file_format = extract_format_from_filename(group_filename)

    # create the folder for group if it doesn't exist
    create_folder(group_path)

    # construct full path to the group fil
    group_path_full = group_path + os.path.sep + group_filename

    # check if the format of each input file match with the group file format
    for path_full in path_full_list:
        input_filename = extract_filename_from_full_path(path_full)
        input_file_format = extract_format_from_filename(input_filename)
        assert (
            input_file_format == group_file_format
        ), f"Format of the group file ({group_file_format}) isn't the same as input file's format ({input_file_format}). Change input file {path_full} or group file {group_path_full}."

    # group input files into a sinlge output file
    with open(group_path_full, "w") as group_file:
        for path_full in path_full_list:
            with open(path_full, "r") as input_file:
                group_file.writelines(input_file)

            group_file.write("\n")

            if delete_input:
                os.remove(path_full)

    return group_path_full


def rescale_density_map(
        density_path_full, 
        rescaled_path_full, 
        box_size=16,
        stderr_file=None
    ):
    """
    Rescales given density map such that it fits the given box size. 
    Achieves this by running relion software commands through a subprocess.

    Params:
    density_path_full - full path to the input denisty map i.e. map to rescale (including its name)
    rescaled_path_full - full path to the output file with rescaled density
    box_size - size of the box (the box to which we rescale)
    stderr_file - file object to log errors from the subprocess if any


    Returns:
    p - object of the Popen class corresponding to the relion's subprocess
    """
  
    # run relion's command through a subpprocess
    p = Popen(
    [
        "relion_image_handler",
        "--i",
        density_path_full,
        "--new_box",
        str(box_size),
        "--o",
        rescaled_path_full
    ],
    stderr=stderr_file,
    )
    
    return p



def log(
    message,
    status="ERROR",
    log_path=os.getcwd() + os.path.sep + "logs",
    log_filename="log.txt",
):
    """
    Creates logs.

    Params:
    message - message to put in the log file
    status - status of the message. Two options: INFO and ERROR
    log_path - path to the directory with log files
    log_filename - name of the log file
    """

    # create directory for log files if it doesn't exist
    # NOTE: commented out folder creation since there were conflicts with Multi Processing
    # create_folder(log_path)

    # full path to the log file (including its name)
    path_full = log_path + os.path.sep + log_filename

    # convert message to string in the case if exception was provided
    message = str(message)

    if not message.endswith("\n"):
        message += "\n"

    # log_message includes utc time of the message, status of the message
    # and the message itself
    log_message = (
        datetime.now(timezone.utc).strftime("%H:%M:%S.%f") + " " + status + ": " + message
    )

    with open(path_full, "a") as f:
        f.write(log_message)


def create_folder(folder_path):
    """
    Creates a folder if it doesn't exist. Usually used for log folders.

    Params:
    folder_path - path to the folder that shold be created
    """

    if not os.path.exists(folder_path):
        os.mkdir(folder_path)