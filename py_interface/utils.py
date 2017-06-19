import os
import logging
import subprocess


log = logging.getLogger('gramtools')


GRAMTOOLS_INSTALL_PATH = os.path.abspath((os.path.join(__file__, '../..')))
gramtools_exec_fpath = os.path.join(GRAMTOOLS_INSTALL_PATH,
                                    'bin',
                                    'gramtools')
prg_build_exec_fpath = os.path.join(GRAMTOOLS_INSTALL_PATH,
                                    'utils',
                                    'vcf_to_linear_prg.pl')
kmers_script_fpath = os.path.join(GRAMTOOLS_INSTALL_PATH,
                                  'utils', 'variantKmers.py')


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
