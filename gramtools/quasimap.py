import os
import time
import json
import logging
import subprocess
import collections

from . import version
from . import common
from . import paths


log = logging.getLogger('gramtools')


def parse_args(common_parser, subparsers):
    parser = subparsers.add_parser('quasimap',
                                   parents=[common_parser])
    parser.add_argument('--gram-directory',
                        help='',
                        type=str)
    parser.add_argument('--reads',
                        help='',
                        type=str)
    parser.add_argument('--output-directory',
                        help='',
                        type=str,
                        required=False)
    parser.add_argument('--kmer-size',
                        help='',
                        type=int,
                        default=15,
                        required=False)


def _execute_command(quasimap_paths, args):
    command = [
        common.gramtools_exec_fpath,
        'quasimap',
        '--gram', quasimap_paths['project'],
        '--reads', quasimap_paths['reads'],
        '--kmer-size', str(args.kmer_size),
        '--run-directory', quasimap_paths['quasimap_run_dirpath'],
    ]

    command_str = ' '.join(command)
    log.debug('Executing command:\n\n%s\n', command_str)

    current_working_directory = os.getcwd()
    process_handle = subprocess.Popen(command_str,
                                      cwd=current_working_directory,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE,
                                      shell=True)

    command_result, entire_stdout = common.handle_process_result(process_handle)

    log.info('Output run directory:\n%s', quasimap_paths['quasimap_run_dirpath'])
    return command_str, command_result, entire_stdout


def _save_report(command_str,
                 command_result,
                 start_time,
                 entire_stdout,
                 quasimap_paths):

    end_time = str(time.time()).split('.')[0]
    _, report_dict = version.report()

    report = collections.OrderedDict([
        ('start_time', start_time),
        ('end_time', end_time),
        ('total_runtime', int(end_time) - int(start_time)),
        ('version_report', report_dict),
        ('command_return_eq_0', command_result),
        ('entire_stdout', entire_stdout),
        ('command_str', command_str),
        ('paths', quasimap_paths),
    ])

    with open(quasimap_paths['run_report'], 'w') as fhandle:
        json.dump(report, fhandle, indent=4)


def run(args):
    log.info('Start process: quasimap')

    start_time = str(time.time()).split('.')[0]
    quasimap_paths = paths.generate_quasimap_paths(args, start_time)
    paths.check_project_file_structure(quasimap_paths)

    results = _execute_command(quasimap_paths, args)
    command_str, command_result, entire_stdout = results
    log.info('End process: quasimap')

    log.debug('Writing run report to run directory')
    _save_report(command_str,
                 command_result,
                 start_time,
                 entire_stdout,
                 quasimap_paths)
