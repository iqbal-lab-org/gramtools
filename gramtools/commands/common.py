## @file
# File path and process result logging related utilities. Are used across commands.
import os
import hashlib
import logging
import subprocess
import gzip
from typing import List, Dict, NamedTuple
from pathlib import Path
from collections import OrderedDict

from Bio import SeqIO

from gramtools import gramtools_lib_fpath

log = logging.getLogger("gramtools")


# Add the dynamically linked libraries at runtime
_old_ld_library_path = (
    os.environ["LD_LIBRARY_PATH"] if "LD_LIBRARY_PATH" in os.environ else ""
)
lib_paths = gramtools_lib_fpath + ":" + _old_ld_library_path


class CommandResult(NamedTuple):
    success: bool
    error_code: int
    stdout: str = ""
    stderr: str = ""


def run_subprocess(command: List) -> CommandResult:
    log.debug("Executing command:\n\n%s\n", " ".join(command))

    current_working_directory = os.getcwd()
    log.debug("Using current working directory:\n%s", current_working_directory)

    process_handle = subprocess.Popen(
        command,
        cwd=current_working_directory,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env={"LD_LIBRARY_PATH": lib_paths},
        universal_newlines=True,
    )

    command_result = handle_process_result(process_handle)
    return command_result


def handle_process_result(process_handle: subprocess.Popen) -> CommandResult:
    """Report process results to logging."""
    captured_stdout = ""
    for line in iter(process_handle.stdout):
        captured_stdout += f"{line}"
        print(line, end="")

    stdout, stderr = process_handle.communicate()

    error_code = process_handle.returncode
    if error_code != 0:
        return CommandResult(False, error_code, captured_stdout, stderr)
    else:
        return CommandResult(True, error_code, captured_stdout, stderr)


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


Seq = str
Chroms = Dict[str, Seq]


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


def write_coordinates_file(chrom_seqs: Chroms, out_fname: str):
    with open(out_fname, "w") as fhandle_out:
        for name, seq in chrom_seqs.items():
            fhandle_out.write(f"{name}\t{len(seq)}\n")
