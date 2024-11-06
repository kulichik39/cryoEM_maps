import os
import sys
import numpy as np
import getopt
from datetime import datetime

# sys.dont_write_bytecode = (
#     True  # disable generating .pyc compiled bytecode files for imported modules
# )

from chimera_utilities import log, save_chimera_density_to_file_with_log
from chimera_utilities import run_chimera_command_with_log as run_com


# create name for the start log file which will include import errors and script's arguments errors if any
# this log file is not related to the model defined by a .pdb file and won't be created if there are no import
# or argument errors
# log file related to the .pdb model will be created later in the code if there are no import
# or argument errors
log_fname_start = datetime.utcnow().strftime("%d_%m_%Y_%H.%M.%S.%f")[:-2] + "_log.txt"


# parse script's arguments
try:
    # i stands for full path to the input .pdb file with model data
    opts, args = getopt.getopt(sys.argv[1:], "i:")

    pdb_path = None  # full path to the model's pdb file (including file's name)
    for opt, arg in opts:
        if opt == "-i":
            pdb_path = arg

    # check if pdb_path is not None after arguments parsing
    assert pdb_path, "Path to the input .pdb file is None after argument's parsing."

    # extract name of the .pdb file
    pdb_fname = pdb_path.split(os.path.sep)[-1]

    # check if the model file indeed ends with .pdb
    assert pdb_fname.endswith(".pdb"), "Name of the model file must end with .pdb ."

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
# that is related to the given .pdb model file
log_fname = (
    datetime.utcnow().strftime("%d_%m_%Y_%H.%M.%S.%f")[:-2]
    + "_"
    + pdb_fname
    + "_log.txt"
)

# run Chimera commands
log("Started Chimera commands.", status="INFO", log_filename=log_fname)

run_com("open " + pdb_path, log_filename=log_fname)  # open model from .pdb file
run_com("molmap #0 1", log_filename=log_fname)  # generate density map

# save density map to file
density_fname = "density_map" + "_" + pdb_fname + ".txt"
save_chimera_density_to_file_with_log(
    density_filename=density_fname, log_filename=log_fname
)

# if all commands
log(
    "Successfully finished Chimera commands! Exiting...",
    status="INFO",
    log_filename=log_fname,
)

run_com("stop now", log_filename=log_fname)
