import os
import sys
import getopt
from datetime import datetime

# sys.dont_write_bytecode = (
#     True  # disable generating .pyc compiled bytecode files for imported modules
# )

from chimera_utilities import (
    log,
    delete_extension_from_filename,
)
from chimera_utilities import run_chimera_command as run_com

# parse script's arguments
opts, args = getopt.getopt(sys.argv[1:], "i:r:o:l:")

molecule_path_full = None  # full path to the molecule file (including its name)
density_path_full = None # full path to the output density file (including its name)
density_resolution = None  # density map resolution (in Angrstrom)
is_log = False # whether we should write logs for Chimera script
log_path = None # path to the log folder (excluding log file name)
log_fname = None # name of the log file
for opt, arg in opts:
    if opt == "-i":
        molecule_path_full = arg
    elif opt == "-r":
        density_resolution = arg
    elif opt == "-o":
        density_path_full = arg
    elif opt == "-l":
        is_log = True
        log_path = arg

# check if the arguments are not None after arguments parsing
assert (
    molecule_path_full
), "Path to the input molecule file is None after argument's parsing."

assert density_resolution, "Density resolution is None after argument's parsing."
assert density_path_full, "Path to the output density file is None after argument's parsing."

if is_log:
    assert log_path, "Path to the log folder is None after argument's parsing."

# extract name of the molecule file
molecule_fname = molecule_path_full.split(os.path.sep)[-1]

# create log file name if required
if is_log:
    log_fname = (
        datetime.utcnow().strftime("%d_%m_%Y_%H.%M.%S.%f")[:-2]
        + "_"
        + delete_extension_from_filename(molecule_fname)
        + "_log.txt"
    )

# run Chimera commands
if is_log:
    log("Started Chimera commands.", status="INFO", log_path=log_path, log_filename=log_fname)

run_com(
    "open " + molecule_path_full, is_log=is_log, log_path=log_path, log_filename=log_fname
)  # open model from the molecule file

# generate density map
volume_id = 0
run_com(
    "molmap #0 {} modelId {}".format(density_resolution, volume_id),
    is_log=is_log,
    log_path=log_path,
    log_filename=log_fname,
)  


# save density map to .mrc file
run_com(
    "volume #{} save {} saveFormat mrc".format(volume_id, density_path_full),
    is_log=is_log,
    log_path=log_path,
    log_filename=log_fname,
)  

# if all commands run successfully, log info message and stop Chimera
if is_log:
    log(
        "Successfully finished Chimera commands! Exiting...",
        status="INFO",
        log_path=log_path,
        log_filename=log_fname,
    )

run_com(
    "stop now",
    is_log=is_log,
    log_path=log_path,
    log_filename=log_fname,
)  