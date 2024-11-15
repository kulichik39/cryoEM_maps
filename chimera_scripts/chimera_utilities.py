import os
import numpy as np
from datetime import datetime
from chimera import runCommand as rc  # use 'rc' as shorthand for runCommand
from chimera import openModels


def log(
    message,
    status="ERROR",
    log_path=os.getcwd() + os.path.sep + "chimera_logs",
    log_filename="log.txt",
):
    """
    Creates logs for python scripts run inside Chimera.

    Params:
    message - message to put in the log file
    status - status of the message. Two options: INFO and ERROR
    log_path - path to the directory with log files
    log_filename - name of the log file
    """

    # create directory for log files if it doesn't exist
    # NOTE: commented out folder creation since there were conflicts with Multi Processing
    # if not os.path.exists(log_path):
    #     os.mkdir(log_path)

    # full path to the log file (including its name)
    path_full = log_path + os.path.sep + log_filename

    # convert message to string in the case if exception was provided
    message = str(message)

    if not message.endswith("\n"):
        message += "\n"

    # log_message includes utc time of the message, status of the message
    # and the message itself
    log_message = (
        datetime.utcnow().strftime("%H:%M:%S.%f") + " " + status + ": " + message
    )

    with open(path_full, "a") as f:
        f.write(log_message)


def run_chimera_command(command):
    """
    Runs chimera's command with runCommand function.

    Params:
    command - string representing a command to chimera
    """

    rc(command)


def run_chimera_command_with_log(
    command, log_path=os.getcwd() + os.path.sep + "chimera_logs", log_filename="log.txt"
):
    """
    Same as run_chimera_command() but with logging options.

    Params:
    command - string representing a command to chimera
    log_path - path to the directory with log files
    log_filename - name of the log file
    """

    try:
        run_chimera_command(command)
        log(
            "Successfully executed chimera command '{}'.".format(command),
            status="INFO",
            log_path=log_path,
            log_filename=log_filename,
        )

    except Exception as e:
        log(
            "Failed to execute chimera command '{}': {}.".format(command, e),
            status="ERROR",
            log_path=log_path,
            log_filename=log_filename,
        )
        raise Exception(e)


def delete_extension_from_filename(filename):
    """
    Deletes extension (symbols after .) from the given filename.

    Params:
    filename - name of the given file

    Returns:
    name of the file without extension
    """

    return ".".join(filename.split(".")[:-1])


def save_chimera_density_to_mrc_file(
    density_path_full,
    volume_id=0
):
    """
    Saves density map computed in Chimera to the specified .mrc file.

    Params:
    density_path_full - full path to the density file (including its name)
    volume_id - id of the Chimera's volume with density data
    """
    
    assert density_path_full.endswith(".mrc"), "mrc filename must end with .mrc!"

    run_chimera_command(
        "volume #{} save {} saveFormat mrc".format(volume_id, density_path_full)
    )


def save_chimera_density_to_mrc_file_with_log(
    density_path_full,
    volume_id=0,
    log_path=os.getcwd() + os.path.sep + "chimera_logs",
    log_filename="log.txt",
):
    """
    Same as save_chimera_density_to_mrc_file() but with logging options.

    Params:
    density_path_full - full path to the density file (including its name)
    volume_id - id of the Chimera's volume with density data
    log_path - path to the directory with log files
    log_filename - name of the log file
    """
    try:
        save_chimera_density_to_mrc_file(
            density_path_full,
            volume_id=volume_id
        )
        log(
            "Successfully saved Chimera density to .mrc file.",
            status="INFO",
            log_path=log_path,
            log_filename=log_filename,
        )
    except Exception as e:
        log(
            "Failed to save Chimera density to .mrc file: {}.".format(e),
            status="ERROR",
            log_path=log_path,
            log_filename=log_filename,
        )
        raise Exception(e)


# The functions save_chimera_density_to_txt_file... save densisty map data to a .txt file.
# They probably won't be used in the code because we usually want to save to the .mrc format,
# which can be achieved by calling save_chimera_density_to_mrc_file. But I keep them just in
# case.
def save_chimera_density_to_txt_file(
    density_path=os.getcwd() + os.path.sep + "density_maps",
    density_filename="density.txt",
):
    """
    Saves density map computed in Chimera to the specified .txt file.

    Params:
    density_path - path to the directory where density maps are stored
    density_fname - name of the file to store the density map in
    """

    assert density_filename.endswith(".txt"), "txt filename must end with .txt!"

    # get density matrix
    _, v = openModels.list()
    data = v.data
    matrix = data.matrix()

    # create directory for density maps if it doesn't exist
    if not os.path.exists(density_path):
        os.mkdir(density_path)

    # full path to the denisty file (including its name)
    path_full = density_path + os.path.sep + density_filename

    with open(path_full, "w") as file:
        file.write("# Array shape:{}\n".format(tuple(int(dim) for dim in matrix.shape)))
        for m_slice in matrix:
            np.savetxt(file, m_slice)
            file.write("# New slice\n")


def save_chimera_density_to_txt_file_with_log(
    density_path=os.getcwd() + os.path.sep + "density_maps",
    density_filename="density.txt",
    log_path=os.getcwd() + os.path.sep + "chimera_logs",
    log_filename="log.txt",
):
    """
    Same as save_chimera_density_to_txt_file() but with logging options.

    Params:
    density_path - path to the directory where density maps are stored
    density_filename - name of the file to store the density map in
    log_path - path to the directory with log files
    log_filename - name of the log file
    """

    try:
        save_chimera_density_to_txt_file(
            density_path=density_path, density_filename=density_filename
        )
        log(
            "Successfully saved Chimera density to .txt file.",
            status="INFO",
            log_path=log_path,
            log_filename=log_filename,
        )
    except Exception as e:
        log(
            "Failed to save Chimera density to .txt file: {}.".format(e),
            status="ERROR",
            log_path=log_path,
            log_filename=log_filename,
        )
        raise Exception(e)
