import os
import hashlib
import logging


log = logging.getLogger('gramtools')

base_install_path = os.path.dirname(os.path.abspath(__file__))
gramtools_exec_fpath = os.path.join(base_install_path,
                                    'bin',
                                    'gram')
prg_build_exec_fpath = os.path.join(base_install_path,
                                    'utils',
                                    'vcf_to_linear_prg.pl')

_old_ld_library_path = os.environ["LD_LIBRARY_PATH"] if "LD_LIBRARY_PATH" in os.environ else ''
lib_paths = (os.path.join(base_install_path, 'lib')
             + ":" + _old_ld_library_path)


def handle_process_result(process_handle):
    """Report process results to logging."""
    if process_handle.stdout is None:
        return True

    uses_stdout = False
    entire_stdout = []
    for line in iter(process_handle.stdout.readline, b''):
        if not uses_stdout:
            log.info('stdout:\n')
            uses_stdout = True
        formatted_line = line.decode('ascii')[:-1]
        print(formatted_line)
        entire_stdout.append(formatted_line)
    if uses_stdout:
        print('')

    stdout, termination_message = process_handle.communicate()
    error_code = process_handle.returncode

    if termination_message:
        log.info('Process termination message:\n%s',
                 termination_message.decode("utf-8"))
        log.info('Process termination code: %s', error_code)

    if stdout:
        log.info('stdout:\n%s', stdout.decode("utf-8"))

    if error_code != 0:
        log.error('Error code != 0')
        return False, entire_stdout
    return True, entire_stdout


def hash_command_paths(command_paths):
    command_hash_paths = {}
    for command, path_component in command_paths.items():
        if isinstance(path_component, list):
            paths = path_component
            command_hash_paths[command] = {p: _file_hash(p) for p in paths if os.path.isfile(p)}
            continue

        path = path_component
        if not os.path.isfile(path):
            continue
        command_hash_paths[command] = _file_hash(path)
    return command_hash_paths


def _file_hash(file_path):
    sha = hashlib.sha256()
    bytes_buffer_size = int(1e+7)
    with open(file_path, 'rb') as f:
        while True:
            data = f.read(bytes_buffer_size)
            if not data:
                break
            sha.update(data)
    return sha.hexdigest()
