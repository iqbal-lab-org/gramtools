import os
import logging
import subprocess


GRAMTOOLS_INSTALL_PATH = '/home/robyn/Documents/gramtools'
gramtools_exec_fpath = os.path.join(GRAMTOOLS_INSTALL_PATH,
                                    'cmake-build-debug',
                                    'gramtools')
prg_build_exec_fpath = os.path.join(GRAMTOOLS_INSTALL_PATH,
                                    'utils',
                                    'vcf_to_linear_prg.pl')
kmers_script_fpath = os.path.join(GRAMTOOLS_INSTALL_PATH,
                                  'utils', 'variantKmers.py')


def setup_logging():
    log = logging.getLogger('gramtools')
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(logging.DEBUG)
    return log

log = setup_logging()


def check_path_exist(paths):
    log.debug('Checking file paths')
    missing_paths = set(path for path in paths
                        if not os.path.isfile(path))

    all_paths_present = not missing_paths
    if all_paths_present:
        return

    log.error("The following paths do not exist:")
    for path in missing_paths:
        log.error(path)
    exit(-1)


def handle_process_result(process_handle):
    """Report process results to logging."""
    stdout, error_message = process_handle.communicate()
    error_code = process_handle.returncode

    if error_message:
        log.info('Process error message:\n%s',
                 error_message.decode("utf-8"))
        log.info('Process error code: %s', error_code)

    if stdout:
        truncate_stdout = stdout.decode("utf-8")
        log.info('stdout:\n%s', truncate_stdout)
    else:
        log.info('stdout: none or piped out')

    if error_code != 0:
        log.error('Error code != 0')
        return False
    return True
