import os
import time
import json
import logging
import subprocess
import collections

from . import common
from . import paths

try:
    from .version import version
except ImportError:
    from .version import fallback_version as version

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
        '--prg', quasimap_paths['prg'],
        '--encoded-prg', quasimap_paths['encoded_prg'],
        '--fm-index', quasimap_paths['fm_index'],
        '--variant-site-mask', quasimap_paths['variant_site_mask'],
        '--allele-mask', quasimap_paths['allele_mask'],
        '--memory-log', quasimap_paths['sdsl_memory_log'],
        '--kmers-prefix-diffs', quasimap_paths['kmer_prefix_diffs'],
        '--kmer-index', quasimap_paths['kmer_index'],
        '--kmer-size', str(args.kmer_size),

        '--reads', quasimap_paths['reads'],
        '--allele-coverages', quasimap_paths['allele_coverage'],
        '--reads-progress', quasimap_paths['reads_progress'],
    ]

    callgrind_command = [
        'valgrind',
        '--tool=callgrind',
        '--callgrind-out-file=' + quasimap_paths['callgrind_out'],
    ]

    callgrind_command = callgrind_command if args.profile else []
    command = callgrind_command + command

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
    commits = version.commit_log.split('*****')[1:]
    commits = '\n'.join(commits)

    end_time = str(time.time()).split('.')[0]

    report = collections.OrderedDict([
        ('start_time', start_time),
        ('end_time', end_time),
        ('total_runtime', int(end_time) - int(start_time)),
        ('current_git_branch', version.current_branch),
        ('command_return_eq_0', command_result),
        ('entire_stdout', entire_stdout),
        ('command_str', command_str),
        ('latest_commit_hash', version.latest_commit),
        ('truncated_commit_log', commits),
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
