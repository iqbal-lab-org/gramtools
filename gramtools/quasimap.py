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

    execute_report = collections.OrderedDict([
        ('command', command_str),
        ('return_value_is_0', command_result),
        ('stdout', entire_stdout),
    ])
    return execute_report


def _save_report(start_time,
                 execute_reports,
                 command_paths,
                 report_file_path):

    end_time = str(time.time()).split('.')[0]
    _, report_dict = version.report()
    current_working_directory = os.getcwd()

    report = collections.OrderedDict([
        ('start_time', start_time),
        ('end_time', end_time),
        ('total_runtime', int(end_time) - int(start_time)),
    ])
    report.update(execute_reports)
    report.update(collections.OrderedDict([
        ('current_working_directory', current_working_directory),
        ('paths', command_paths),
        ('version_report', report_dict),
    ]))

    with open(report_file_path, 'w') as fhandle:
        json.dump(report, fhandle, indent=4)


def run(args):
    log.info('Start process: quasimap')

    start_time = str(time.time()).split('.')[0]
    quasimap_paths = paths.generate_quasimap_paths(args, start_time)
    paths.check_project_file_structure(quasimap_paths)

    gramtools_cpp_report = _execute_command(quasimap_paths, args)

    log.debug('Writing run report to run directory')
    execute_reports = collections.OrderedDict([
        ('gramtools_cpp_quasimap', gramtools_cpp_report),
    ])
    report_file_path = quasimap_paths['run_report']
    _save_report(start_time,
                 execute_reports,
                 quasimap_paths,
                 report_file_path)
    log.info('End process: quasimap')
