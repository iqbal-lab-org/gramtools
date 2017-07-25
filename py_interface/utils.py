import os
import logging


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
