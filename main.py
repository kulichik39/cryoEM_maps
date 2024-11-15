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
    log
)

from internal_flexibility import (
    generate_conformers,
)

from missing_parts import random_delete_atoms_from_pdb_file

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
        use_small_ring_torsions=True, 
        prune_rms_tresh=1.0,
        random_seed=0xF00D,
        delete_prob=0.2,
        density_resolution=3.5,
        box_size=16
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
    """

    # create the folder for temporary files if it doesn't exist
    if not os.path.exists(temporary_path): 
        os.mkdir(temporary_path)

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
    p = compute_density_map_in_chimera(
        delatoms_conformer_path, density_path_full, density_resolution=density_resolution
    )
    _, stderr = p.communicate()
    # if the Chimera's subprocess finished with errors, raise RuntimeError
    if stderr:
        raise RuntimeError(f"Chimera's subprocess finished with errors: {stderr}")
    

    # rescale density map using Relion software
    p = rescale_density_map(density_path_full, output_density_path_full, box_size=box_size)
    _, stderr = p.communicate()

    # if the Relion's subprocess finished with errors, raise RuntimeError
    if stderr:
        raise RuntimeError(f"Relion's subprocess finished with errors: {stderr}")


def main(
        molecule_name,
        db_path,
        temporary_path=os.getcwd() + os.path.sep + "temp",
        remove_Hs=True, 
        n_confs=5, 
        add_Hs=False, 
        use_small_ring_torsions=True, 
        prune_rms_tresh=1.0,
        random_seed=0xF00D,
        delete_prob=0.2,
        density_resolution=3.5,
        box_size=16,
        write_log=True,
        log_filename="log.txt",
        log_path=os.getcwd() + os.path.sep + "log",
        clear_temporary=True
        ):
    """
    The main function for generating low resolution density maps for molecule ligands from the database.

    Params:
    molecule_name - name of the molecule
    db_path - path to the database
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
    write_log - whether logs should be written
    log_filename - name of the log file
    log_path - path to the folder where log file should be written (excluding its name)
    clear_temporary - whether to clear folder with temporary files once the functions exits

    Returns:
    output_density_path_full - full path to the output density file (including its name)
    """

    try:
        # find corresponding ligand path in the database
        input_ligand_path_full = find_ligand_file_in_db(molecule_name, db_path)

        # construct full path to the output density file
        output_density_path_full = (
            db_path + 
            os.path.sep + 
            molecule_name + 
            os.path.sep + 
            f"{molecule_name}_ligand_scaled_boxed_{box_size}A_low_resolution_forward_model.mrc"
        )
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
            box_size=box_size
        )

        if write_log:
            log(f"Successfully computed density map for {molecule_name}. Check the map in: {output_density_path_full}.", status="INFO", log_filename=log_filename, log_path=log_path)

    except Exception as e:
        if write_log:
            log(f"Failed to compute density map for {molecule_name}: {e}", status="ERROR", log_filename=log_filename, log_path=log_path)

    finally:
        if clear_temporary:
            shutil.rmtree(temporary_path)

    return output_density_path_full



if __name__ == "__main__":

    # path to the database with molecule data
    db_path = os.path.sep + os.path.sep.join(["mnt", "cephfs", "projects", "2023110101_Ligand_fitting_to_EM_maps", "PDBbind", "PDBBind_Zenodo_6408497"])

    # load molecule names
    molecule_names_csv = db_path + os.path.sep + "PDB_IDs_with_rdkit_length_less_than_16A_succ_gnina.csv" # path to the .csv file with molecule names
    molecule_names = read_molecule_names_from_file(molecule_names_csv)


    # create log folders
    log_path = os.getcwd() + os.path.sep + "log" # logs for main.py
    if not os.path.exists(log_path):
        os.mkdir(log_path)

    chimera_log_path = os.getcwd() + os.path.sep + "chimera_logs" # logs for chimera scripts
    if not os.path.exists(chimera_log_path):
            os.mkdir(chimera_log_path)

    # run computations in several Processes to speed up
    procs = []
    for molecule_name in molecule_names:
        temporary_path = os.getcwd() + os.path.sep + f"{molecule_name}_temp" # path to the folder where temporary files will be stored (e.g. conformers)
        kwargs = {"temporary_path": temporary_path, "clear_temporary": True, "log_path": log_path}
        proc = Process(target=main, args=(molecule_name, db_path), kwargs=kwargs)
        procs.append(proc)
        proc.start()

    for proc in procs:
        proc.join()



