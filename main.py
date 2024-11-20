import os
import numpy as np
import shutil
from multiprocessing import Process
from utilities import (
    read_molecule,
    delete_extension_from_filename,
    save_all_conformers_to_pdb,
    extract_filename_from_full_path,
    group_conformers_to_single_file,
    compute_density_map_in_chimera,
    rescale_density_map,
    log,
    create_folder
)

from internal_flexibility import (
    generate_conformers,
)

from missing_parts import random_delete_atoms_from_pdb_file
from datetime import datetime, timezone

def read_molecule_names_from_file(path_full, skip_header=1, delimiter=','):
    """
    Reads molecules' names from the given file and stores them in numpy array.

    Params:
    path_full - full path to the input file
    skip_header - how many header rows to skip when read the file
    delimiter - delimiter used in the input file (e.g. ',' for .csv)

    Returns:
    molecule_names - numpy array with molecule names    
    """

    molecule_names = np.genfromtxt(path_full, skip_header=skip_header, delimiter=delimiter, dtype=np.str_)

    assert len(molecule_names.shape) == 1, f"Array with molecule names should be 1D, but has the shape: {molecule_names.shape}"

    return molecule_names

def find_ligand_file_in_db(molecule_name, db_path):
    """
    Finds .pdb or .sdf file for the given molecule's ligand in the database.
    At first, tries to look for a .pdb file and if found returns full path to the file (icnluding its name).
    Then, tries to look for a .sdf file and if found returns full path to the file (icnluding its name).
    If none of the files was found, raises RuntimeError.

    Params:
    molecule_name - name of the molecule
    db_path - path to the database with molecule files

    Returns:
    full path to the ligand file (icnlduing its name)
    """

    # path to the folder with the given molecule data
    molecule_folder = db_path + os.path.sep + molecule_name

    # try to find .pdb file
    ligand_filename = molecule_name + "_ligand.pdb"
    if os.path.isfile(molecule_folder + os.path.sep + ligand_filename):
        return molecule_folder + os.path.sep + ligand_filename

    # if .pdb wasn't found, try to find .sdf file
    ligand_filename = molecule_name + "_ligand.sdf"
    if os.path.isfile(molecule_folder + os.path.sep + ligand_filename):
        return molecule_folder + os.path.sep + ligand_filename
    
    raise RuntimeError(f"Ligand file for the {molecule_name} molecule wasn't found!")



def generate_low_resolution_density(
        input_molecule_path_full,
        output_density_path_full, 
        temporary_path=os.getcwd() + os.path.sep + "temp",
        remove_Hs=True, 
        n_confs=5, 
        add_Hs=False, 
        use_small_ring_torsions=False, 
        prune_rms_tresh=1.0,
        random_seed=0xF00D,
        delete_prob=0.2,
        density_resolution=3.5,
        box_size=16,
        is_log_chimera=False,
        log_chimera_path=os.getcwd() + os.path.sep + "chimera_logs"
        ):
    """
    Creates low resolution density map for the given molecule.

    Params:
    input_molecule_path_full - full path to the input file with molecule data (including its name)
    output_density_path_full - full path to the output file with density data (including its name)
    temporary_path - path to the folder where temporary files (e.g conformers) will be stored
    remove_Hs - whether to remove hydrogen atoms when reading  the input molecule file
    n_confs - number of conformers to generate
    add_Hs - whether to add hydrogen atoms to the given molecule before conformers generation
    use_small_ring_torsions - whether to include additional small ring torsion potentials
    prune_rms_tresh - threshold value for RMSD pruning of the conformers
    random_seed - seed for conformers generation algorithm (for reproducibility)
    delete_prob - probability of deleting an atom from the molecule file
    density_resolution - desired resolution of the density map (in Angstrom)
    box_size - size of the box for density rescaling
    is_log_chimera - whether logs for Chimera script should be written
    log_chimera_path - path to the folder where Chimera log files will be written
    """

    # create the folder for temporary files if it doesn't exist
    create_folder(temporary_path)

    # read molecule data using RDKit
    mol = read_molecule(input_molecule_path_full, remove_Hs=remove_Hs)

    # generate conformers
    mol = generate_conformers(
        mol,
        n_confs=n_confs,
        add_Hs=add_Hs,
        use_small_ring_torsions=use_small_ring_torsions,
        prune_rms_tresh=prune_rms_tresh,
        random_seed=random_seed
    )

    # save generated conformers to files
    molecule_filename = extract_filename_from_full_path(input_molecule_path_full)
    base_conformer_filename = delete_extension_from_filename(molecule_filename) + ".pdb"
    conformers_folder = temporary_path + os.path.sep + "raw_conformers_data"
    conformer_path_list = save_all_conformers_to_pdb(
        mol, base_conformer_filename, pdb_path=conformers_folder 
    )

    # group conformers' files into a single file to generate an "average" map
    group_filename = (
        "group_conformers_" + base_conformer_filename
    )
    group_conformers_path = group_conformers_to_single_file(
        conformer_path_list,
        group_filename,
        group_path=conformers_folder,
        delete_input=False,
    )

    # randomly delete atoms from the generated conformers
    delatoms_conformer_folder = temporary_path + os.path.sep  + "delatoms_conformers_data"
    delatoms_conformer_path = random_delete_atoms_from_pdb_file(
        group_conformers_path,
        delatoms_molecule_path=delatoms_conformer_folder,
        delete_prob=delete_prob,
    )


    # compute "average" density map for the grouped conformers using Chimera software
    density_filename = (
            delete_extension_from_filename(molecule_filename)
            + "_low_resolution_forward_model"
            + ".mrc"
        )
    density_path_full = temporary_path + os.path.sep + density_filename

    # full path to the file to store Chimera subprocess errors if any
    stderr_file_chim_path = (
        os.getcwd() + 
        os.path.sep + 
        f"temp_chim_{delete_extension_from_filename(extract_filename_from_full_path(input_molecule_path_full))}.err"
    )
    stderr_file_chim = open(stderr_file_chim_path, "w") # file to store Chimera errors 
    p_chim = compute_density_map_in_chimera(
        delatoms_conformer_path, 
        density_path_full, 
        is_log=is_log_chimera, 
        log_path=log_chimera_path, 
        density_resolution=density_resolution,
        stderr_file=stderr_file_chim
    )
    return_code_chim = p_chim.wait()
    stderr_file_chim.close()

    # if the Chimera's subprocess finished with errors, raise RuntimeError
    stderr_file_chim = open(stderr_file_chim_path, "r")
    err_chim = stderr_file_chim.read()
    stderr_file_chim.close()
    os.remove(stderr_file_chim_path)
    # NOTE: here we skip errors that contain "Cannot find consistent set of bond" cause even with it
    # Chimera generates code density maps
    if return_code_chim != 0 or (err_chim and "Cannot find consistent set of bond" not in err_chim):
        raise RuntimeError("Chimera's subprocess finished with errors: " + err_chim)
    

    # rescale density map using Relion software
    # full path to the file to store Relion subprocess errors if any
    stderr_file_relion_path = (
        os.getcwd() + 
        os.path.sep + 
        f"temp_relion_{delete_extension_from_filename(extract_filename_from_full_path(input_molecule_path_full))}.err"
    )
    stderr_file_relion = open(stderr_file_relion_path, "w") # file to store Chimera errors 
    p_relion = rescale_density_map(
        density_path_full, 
        output_density_path_full, 
        box_size=box_size,
        stderr_file=stderr_file_relion
    )
    # _, stderr_relion = p_relion.communicate()
    return_code_relion = p_relion.wait()
    stderr_file_relion.close()

    # if the Relion's subprocess finished with errors, raise RuntimeError
    stderr_file_relion = open(stderr_file_relion_path, "r")
    err_relion = stderr_file_relion.read()
    stderr_file_relion.close()
    os.remove(stderr_file_relion_path)
    if return_code_relion != 0 or err_relion:
        # raise RuntimeError(f"Relion's subprocess finished with errors: {stderr_relion}")
        raise RuntimeError(f"Relion's subprocess finished with errors: " + err_relion)


def main(
        molecule_name,
        db_path,
        output_density_path_full,
        temporary_path=os.getcwd() + os.path.sep + "temp",
        remove_Hs=True, 
        n_confs=5, 
        add_Hs=False, 
        use_small_ring_torsions=False, 
        prune_rms_tresh=1.0,
        random_seed=0xF00D,
        delete_prob=0.2,
        density_resolution=3.5,
        box_size=16,
        clear_temporary=True,
        is_log_main=True,
        log_main_path=os.getcwd() + os.path.sep + "log",
        log_main_filename="log.txt",
        is_log_chimera=False,
        log_chimera_path=os.getcwd() + os.path.sep + "chimera_logs"
        ):
    """
    The main function for generating low resolution density maps for molecule ligands from the database.

    Params:
    molecule_name - name of the molecule
    db_path - path to the database
    output_density_path_full - full path to the output density file (including its name)
    temporary_path - path to the folder where temporary files (e.g conformers) will be stored
    remove_Hs - whether to remove hydrogen atoms when reading  the input molecule file
    n_confs - number of conformers to generate
    add_Hs - whether to add hydrogen atoms to the given molecule before conformers generation
    use_small_ring_torsions - whether to include additional small ring torsion potentials
    prune_rms_tresh - threshold value for RMSD pruning of the conformers
    random_seed - seed for conformers generation algorithm (for reproducibility)
    delete_prob - probability of deleting an atom from the molecule file
    density_resolution - desired resolution of the density map (in Angstrom)
    box_size - size of the box for density rescaling
    clear_temporary - whether to clear folder with temporary files once the functions exits
    is_log_main - whether logs for the main script should be written
    log_main_path - path to the folder where main log file should be written 
    log_main_filename - name of the main log file
    is_log_chimera - whether logs for Chimera script should be written
    log_chimera_path - path to the folder where Chimera log file should be written 
    """

    try:
        # find corresponding ligand path in the database
        input_ligand_path_full = find_ligand_file_in_db(molecule_name, db_path)

        generate_low_resolution_density(
            input_ligand_path_full,
            output_density_path_full,
            temporary_path=temporary_path,
            remove_Hs=remove_Hs,
            n_confs=n_confs,
            add_Hs=add_Hs,
            use_small_ring_torsions=use_small_ring_torsions,
            prune_rms_tresh=prune_rms_tresh,
            random_seed=random_seed,
            delete_prob=delete_prob,
            density_resolution=density_resolution,
            box_size=box_size,
            is_log_chimera=is_log_chimera,
            log_chimera_path=log_chimera_path
        )

        if is_log_main:
            log(f"Successfully computed density map for {molecule_name}. Check the map in: {output_density_path_full}.", status="INFO", 
                log_filename=log_main_filename, log_path=log_main_path)

    except Exception as e:
        if is_log_main:
            log(f"Failed to compute density map for {molecule_name}: {e}", status="ERROR", 
                log_filename=log_main_filename, log_path=log_main_path)

    finally:
        if clear_temporary:
            shutil.rmtree(temporary_path)




if __name__ == "__main__":

    # path to the database with molecule data
    db_path = os.path.sep + os.path.sep.join(["mnt", "cephfs", "projects", "2023110101_Ligand_fitting_to_EM_maps", "PDBbind", "PDBBind_Zenodo_6408497"])

    # load molecule names
    molecule_names_csv = db_path + os.path.sep + "PDB_IDs_with_rdkit_length_less_than_16A_succ_gnina.csv" # path to the .csv file with molecule names
    molecule_names = read_molecule_names_from_file(molecule_names_csv)


    # create log folders. NOTE: important to create if they don't exist, cause the code won't work otherwise!
    is_log_main = True
    is_log_chimera = True
    log_main_path = os.getcwd() + os.path.sep + "log" # logs for main.py
    log_main_filename = (
        datetime.now(timezone.utc).strftime("%d_%m_%Y_%H.%M.%S.%f")[:-2]
        + "_log.txt"
    )
    log_chimera_path = os.getcwd() + os.path.sep + "chimera_logs" # logs for chimera scripts
    if is_log_main:
        create_folder(log_main_path)
    if is_log_chimera:
        create_folder(log_chimera_path)

    # run computations in several Processes to speed up
    n_proc = 94 # number of processes 
    delete_prob = 0.2
    box_size = 16

    i = 0

    # molecule_names = []
    # with open(os.getcwd() + os.path.sep + "log" + os.path.sep + "bad_mols.txt", "r") as f:
    #     for line in f:
    #         molecule_names.append(line.strip())

    # print(molecule_names)


    # molecule_names = molecule_names[10: 1000]
    
    while i < len(molecule_names):
        if i + n_proc <= len(molecule_names):
            last_id = i + n_proc + 1 
        else:
            last_id = len(molecule_names)

        procs = []

        for j in range(i, last_id):
            molecule_name = molecule_names[j]
            output_density_path_full = (
                db_path + 
                os.path.sep + 
                molecule_name + 
                os.path.sep + 
                f"{molecule_name}_ligand_scaled_boxed_{box_size}A_delprob{delete_prob}_low_resolution_forward_model.mrc"
            )
            temporary_path = os.getcwd() + os.path.sep + f"{molecule_name}_temp" # path to the folder where temporary files will be stored (e.g. conformers)
            kwargs = {
                "temporary_path": temporary_path, 
                "clear_temporary": True, 
                "delete_prob": delete_prob, 
                "box_size": box_size,
                "is_log_main": is_log_main,
                "log_main_path": log_main_path,
                "is_log_chimera": is_log_chimera,
                "log_chimera_path": log_chimera_path,
                "log_main_filename": log_main_filename,
                "use_small_ring_torsions": False
            }
            proc = Process(target=main, args=(molecule_name, db_path, output_density_path_full), kwargs=kwargs)
            procs.append(proc)
            proc.start()

        for proc in procs:
            proc.join()

        i = last_id



