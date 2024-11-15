import os
import sys
import getopt
from datetime import datetime

# sys.dont_write_bytecode = (
#     True  # disable generating .pyc compiled bytecode files for imported modules
# )

from chimera_utilities import (
    log,
    save_chimera_density_to_mrc_file_with_log,
    delete_extension_from_filename,
)
from chimera_utilities import run_chimera_command_with_log as run_com


# create name for the start log file which will include import errors and script's arguments
# errors if any
# this log file is not related to the model defined by the input molecule file and
# won't be created if there are no import or argument errors
# log file related to the molecule model will be created later in the code if there are no
# import or argument errors
log_fname_start = datetime.utcnow().strftime("%d_%m_%Y_%H.%M.%S.%f")[:-2] + "_log.txt"


# parse script's arguments
try:
    # i stands for the full path to the input file with molecule data
    opts, args = getopt.getopt(sys.argv[1:], "i:r:o:")

    molecule_path_full = None  # full path to the molecule file (including its name)
    density_path_full = None # full path to the output density file (including its name)
    density_resolution = None  # density map resolution (in Angrstrom)
    for opt, arg in opts:
        if opt == "-i":
            molecule_path_full = arg
        elif opt == "-r":
            density_resolution = arg
        elif opt == "-o":
            density_path_full = arg

    # check if the arguments are not None after arguments parsing
    assert (
        molecule_path_full
    ), "Path to the input molecule file is None after argument's parsing."

    assert density_resolution, "Density resolution is None after argument's parsing."
    assert density_path_full, "Path to the output density file is None after argument's parsing."

    # extract name of the molecule file
    molecule_fname = molecule_path_full.split(os.path.sep)[-1]

except getopt.GetoptError as e:
    log(
        "Failed to parse script's arguments: {}.".format(e),
        status="ERROR",
        log_filename=log_fname_start,
    )
    raise getopt.GetoptError(e)

except AssertionError as e:
    log(
        e,
        status="ERROR",
        log_filename=log_fname_start,
    )
    raise AssertionError(e)


# If all the exceptions above successfully passed, create a new log file
# that is related to the given model file
log_fname = (
    datetime.utcnow().strftime("%d_%m_%Y_%H.%M.%S.%f")[:-2]
    + "_"
    + delete_extension_from_filename(molecule_fname)
    + "_log.txt"
)

# run Chimera commands
log("Started Chimera commands.", status="INFO", log_filename=log_fname)

run_com(
    "open " + molecule_path_full, log_filename=log_fname
)  # open model from the molecule file

volume_id = 0
run_com(
    "molmap #0 {} modelId {}".format(density_resolution, volume_id),
    log_filename=log_fname,
)  # generate density map

# save density map to file
save_chimera_density_to_mrc_file_with_log(
    density_path_full,
    volume_id=volume_id,
    log_filename=log_fname,
)

# if all commands run successfully, log info message and stop Chimera
log(
    "Successfully finished Chimera commands! Exiting...",
    status="INFO",
    log_filename=log_fname,
)

run_com("stop now", log_filename=log_fname)
