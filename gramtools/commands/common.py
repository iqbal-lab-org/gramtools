## @file
# File path and process result logging related utilities. Are used across commands.
import os
import hashlib
import logging
import subprocess
from typing import List
from pathlib import Path
from collections import OrderedDict

from Bio import SeqIO

import gramtools

log = logging.getLogger("gramtools")

# Find executable locations
base_install_path = Path(gramtools.__file__).resolve().parent
gramtools_exec_fpath = str(base_install_path / "bin" / "gram")

# Add the dynamically linked libraries at runtime
_old_ld_library_path = (
    os.environ["LD_LIBRARY_PATH"] if "LD_LIBRARY_PATH" in os.environ else ""
)
lib_paths = str(base_install_path / "lib") + ":" + _old_ld_library_path


def run_subprocess(command: List):
    command_str = " ".join(command)
    log.debug("Executing command:\n\n%s\n", command_str)

    current_working_directory = os.getcwd()
    log.debug("Using current working directory:\n%s", current_working_directory)

    process_handle = subprocess.Popen(
        command_str,
        cwd=current_working_directory,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True,
        env={"LD_LIBRARY_PATH": lib_paths},
    )

    command_result, entire_stdout = handle_process_result(process_handle)
    return command_result, entire_stdout


def handle_process_result(process_handle):
    """Report process results to logging."""
    if process_handle.stdout is None:
        return True

    entire_stdout = []
    for line in iter(process_handle.stdout.readline, b""):
        formatted_line = line.decode("ascii")[:-1]
        entire_stdout.append(formatted_line)

    stdout, termination_message = process_handle.communicate()
    error_code = process_handle.returncode

    if termination_message:
        log.info(
            "Process termination message:\n%s", termination_message.decode("utf-8")
        )
        log.info("Process termination code: %s", error_code)

    if error_code != 0:
        log.error("Error code != 0")
        log.error("stdout:\n")
        print("\n".join(entire_stdout))
        return False, entire_stdout
    else:
        log.info("stdout:\n")
        print("\n".join(entire_stdout))
        return True, entire_stdout


def hash_command_paths(command_paths):
    command_hash_paths = {}
    for command, path_component in command_paths.items():
        if isinstance(path_component, list):
            paths = path_component
            command_hash_paths[command] = {
                str(p): _file_hash(p) for p in paths if p.is_file()
            }
            continue

        if path_component.is_file():
            command_hash_paths[command] = _file_hash(path_component)
    return command_hash_paths


def _file_hash(file_path):
    sha = hashlib.sha256()
    bytes_buffer_size = int(1e7)
    with open(file_path, "rb") as f:
        while True:
            data = f.read(bytes_buffer_size)
            if not data:
                break
            sha.update(data)
    return sha.hexdigest()


def load_fasta(reference_file, sizes_only=False):
    ref_records = OrderedDict()
    # SeqIO: the id of the record is everything between ">" and the first space.
    for seq_record in SeqIO.parse(reference_file, "fasta"):
        ref_records[seq_record.id] = str(seq_record.seq)
        if sizes_only:
            ref_records[seq_record.id] = len(ref_records[seq_record.id])
    return ref_records
