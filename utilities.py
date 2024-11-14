import os
import numpy as np
import mrcfile
from subprocess import PIPE, Popen
from rdkit import Chem


def run_script_inside_chimera(
    script_name,
    script_path=os.getcwd() + os.path.sep + "chimera_scripts",
    option_names=(),
    option_values=(),
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
        stdout=PIPE,
        stderr=PIPE,
    )

    return p


def compute_density_map_in_chimera(
    molecule_path_full,
    density_resolution=1.0,
    script_path=os.getcwd() + os.path.sep + "chimera_scripts",
):
    """
    Computes density map for the given molecule/conformer using Chimera.
    To achieve this, runs python script "chimera_density_map.py" inside Chimera through a
    subprocess.

    Params:
    molecule_path_full - full path to the input molecule file (including its name)
    density_resolution - desired resolution of the map (in Angstrom)
    script_path - path to the folder with python script (excluding its name)

    Returns:
    p - object of the Popen class corresponding to the Chimera's subprocess
    """

    # name of the python script to run inside Chimera
    scipt_name = "chimera_density_map.py"

    # options for the python script
    option_names = ("-i", "-r")
    option_values = (molecule_path_full, density_resolution)

    # returns object corresponding to the subprocess
    p = run_script_inside_chimera(
        scipt_name,
        script_path=script_path,
        option_names=option_names,
        option_values=option_values,
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
#     if not os.path.exists(molecule_path):
#         os.mkdir(molecule_path)

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


def extract_filename_from_full_path(full_path):
    """
    Extracts filename from the full (including the name) path to file.

    Params:
    full_path - full path to the file

    Returns:
    filename extracted from the full path
    """

    return full_path.split(os.path.sep)[-1]


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
    sdf_filename,
    sdf_path=os.getcwd() + os.path.sep + "raw_molecule_data",
    n_mols=1,
    remove_Hs=True,
):
    """
    Reads molecules data from the given .sdf file using RDKit library.

    Params:
    sdf_filename - name of the input sdf file
    sdf_path - path to the folder where the input sdf is stored (excluding its name)
    n_mols - expected number of molecules to read from the file
    remove_Hs - whether to remove hydrogen atoms from the molecule when reading

    Returns:
    mols - list with RDKit molecule objects obtained from the given file
    """

    assert sdf_filename.endswith(".sdf"), "sdf filename must end with .sdf!"

    # construct full path to the given sdf file (including its name)
    sdf_path_full = sdf_path + os.path.sep + sdf_filename

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
    pdb_filename,
    pdb_path=os.getcwd() + os.path.sep + "raw_molecule_data",
    remove_Hs=True,
):
    """
    Reads one molecule data from the given .pdb file using RDKit library.
    NOTE: even if you have several molecules in the given file, RDKit will
    read only the first one.

    Params:
    pdb_filename - name of the input pdb file
    pdb_path - path to the folder where the input pdb is stored (excluding its name)
    remove_Hs - whether to remove hydrogen atoms from the molecule when reading

    Returns:
    mol - RDKit molecule object obtained from the given file
    """

    assert pdb_filename.endswith(".pdb"), "pdb filename must end with .pdb!"

    # construct full path to the given pdb file (including its name)
    pdb_path_full = pdb_path + os.path.sep + pdb_filename

    mol = Chem.MolFromPDBFile(pdb_path_full, removeHs=remove_Hs)

    return mol


def read_molecule(
    filename, file_path=os.getcwd() + os.path.sep + "raw_molecule_data", remove_Hs=True
):
    """
    Reads one molecule data from the given file and converts it to an RDKit molecule object.

    Params:
    filename - name of the input file
    file_path - path to the folder where the input file is stored (excluding its name)
    remove_Hs - whether to remove hydrogen atoms from the molecule when reading

    Returns:
    mol - RDKit molecule object obtained from the given file
    """

    # extract format of the input file
    file_format = extract_format_from_filename(filename)

    match file_format:
        case "sdf":
            mols = read_molecules_from_sdf(
                filename, sdf_path=file_path, n_mols=1, remove_Hs=remove_Hs
            )
            return mols[0]

        case "pdb":
            mol = read_molecule_from_pdb(
                filename, pdb_path=file_path, remove_Hs=remove_Hs
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

    if not os.path.exists(pdb_path):  # create the folder for saving if it doesn't exist
        os.mkdir(pdb_path)

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
    delte_input - whether to delete input conformers' files

    Returns:
    group_path_full - full path to the file where conformers are grouped (including its name)
    """

    # exctract format of the group file
    group_file_format = extract_format_from_filename(group_filename)

    if not os.path.exists(
        group_path
    ):  # create the folder for group if it doesn't exist
        os.mkdir(group_path)

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
