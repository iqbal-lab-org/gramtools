## @file
# File path and process result logging related utilities. Are used across commands.
import os
import hashlib
import logging
import subprocess
import gzip
from typing import List, NamedTuple
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


class CommandResult(NamedTuple):
    success: bool
    error_code: int
    stdout: str = ""
    stderr: str = ""


def run_subprocess(command: List) -> CommandResult:
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

    command_result = handle_process_result(process_handle)
    return command_result


def handle_process_result(process_handle: subprocess.Popen) -> CommandResult:
    """Report process results to logging."""
    stdout, stderr = process_handle.communicate()
    stdout = stdout.decode()
    stderr = stderr.decode()

    error_code = process_handle.returncode
    if error_code != 0:
        return CommandResult(False, error_code, stdout, stderr)
    else:
        log.info(f"stdout:\n{stdout}")
        return CommandResult(True, error_code, stdout, stderr)


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


def load_fasta(reference_file: Path, sizes_only=False):
    if str(reference_file).endswith(".gz"):
        ref_fhandle = gzip.open(str(reference_file), "rt")
    else:
        ref_fhandle = open(str(reference_file))

    ref_records = OrderedDict()
    # SeqIO: the id of the record is everything between ">" and the first space.
    for seq_record in SeqIO.parse(ref_fhandle, "fasta"):
        ref_records[seq_record.id] = str(seq_record.seq)
        if sizes_only:
            ref_records[seq_record.id] = len(ref_records[seq_record.id])
    ref_fhandle.close()
    return ref_records
