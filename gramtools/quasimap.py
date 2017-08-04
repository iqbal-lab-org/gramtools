import os
import time
import json
import logging
import subprocess
import collections

from . import common

try:
    from .version import version
except ImportError:
    from .version import fallback_version as version

log = logging.getLogger('gramtools')


def parse_args(common_parser, subparsers):
    parser = subparsers.add_parser('quasimap',
                                   parents=[common_parser])
    parser.add_argument('--gram-files', help='',
                        type=str)
    parser.add_argument('--fastaq', help='',
                        type=str)
    parser.add_argument('--kmer-size', help='',
                        type=int)


def _get_project_dirpath(prg_fpath):
    prg_fname = os.path.basename(prg_fpath)
    if prg_fname.endswith('.prg'):
        project_dir = prg_fname[:-len('.prg')]
    else:
        project_dir = prg_fname
    project_dir = project_dir.split('.')[0]

    current_working_directory = os.getcwd()
    project_dirpath = os.path.join(current_working_directory, project_dir)
    return project_dirpath


def _get_run_dirpath(output_dirpath, project_dirpath,
                     ksize, start_time):
    project = os.path.basename(project_dirpath)

    template = '{time}_{project}_ksize{ksize}'
    run_dir = template.format(time=start_time, project=project,
                              ksize=ksize)

    run_dirpath = os.path.join(output_dirpath, run_dir)
    return run_dirpath


def _get_paths(args, start_time):
    project = os.path.abspath(args.gram_files)
    output_dirpath = project + '_output'

    run_dirpath = _get_run_dirpath(output_dirpath, project,
                                   args.kmer_size, start_time)

    project_root = {
        'project': project,
        'prg': os.path.join(project, 'prg'),
        'sites_mask': os.path.join(project, 'sites_mask'),
        'allele_mask': os.path.join(project, 'allele_mask'),
    }

    kmer_paths = {
        'kmer': os.path.join(project, 'kmer'),
        'kmer_file': os.path.join(project, 'kmer',
                                  'ksize_' + str(args.kmer_size)),
    }

    cache_paths = {
        'cache': os.path.join(project, 'cache'),
        'int_encoded_prg': os.path.join(
            project, 'cache', 'int_encoded_prg'),
        'fm_index': os.path.join(project, 'cache', 'fm_index'),
        'kmer_suffix_array': os.path.join(project, 'cache',
                                          'kmer_suffix_array'),
    }

    output_paths = {
        'output': output_dirpath,
        'run': run_dirpath,
        'run_report': os.path.join(run_dirpath, 'report.json'),
        'callgrind_out': os.path.join(run_dirpath, 'callgrind.out'),
        'fm_index_memory_log': os.path.join(
            run_dirpath, 'fm_index_memory_log'),
        'allele_coverage': os.path.join(run_dirpath, 'allele_coverage'),
        'reads': os.path.join(run_dirpath, 'reads'),
    }

    other_paths = {
        'fastaq': args.fastaq,
    }

    paths = {}
    paths.update(project_root)
    paths.update(kmer_paths)
    paths.update(cache_paths)
    paths.update(output_paths)
    paths.update(other_paths)
    return paths


def _setup_file_structure(paths):
    dirs = [
        paths['project'],
        paths['cache'],
        paths['kmer'],
        paths['kmer_suffix_array'],
        paths['output'],
        paths['run'],
    ]
    for dirpath in dirs:
        if os.path.isdir(dirpath):
            continue
        log.debug('Creating directory:\n%s', dirpath)
        os.mkdir(dirpath)


def _execute_command(paths, args):
    command = [
        common.gramtools_exec_fpath,
        '--prg', paths['prg'],
        '--csa', paths['fm_index'],
        '--ps', paths['sites_mask'],
        '--pa', paths['allele_mask'],
        '--co', paths['allele_coverage'],
        '--ro', paths['reads'],
        '--po', paths['int_encoded_prg'],
        '--log', paths['fm_index_memory_log'],
        '--kfile', paths['kmer_file'],
        '--input', paths['fastaq'],
        '--ksize', str(args.kmer_size),
    ]

    callgrind_command = [
        'valgrind',
        '--tool=callgrind',
        '--callgrind-out-file=' + paths['callgrind_out'],
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

    log.info('Output run directory:\n%s', paths['run'])
    return command_str, command_result, entire_stdout


def _save_report(command_str, command_result, start_time, entire_stdout, paths):
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
        ('paths', paths),
    ])

    with open(paths['run_report'], 'w') as fhandle:
        json.dump(report, fhandle, indent=4)


def run(args):
    log.info('Start process: quasimap')

    start_time = str(time.time()).split('.')[0]
    paths = _get_paths(args, start_time)
    _setup_file_structure(paths)

    command_str, command_result, entire_stdout = _execute_command(paths, args)
    log.info('End process: quasimap')

    log.debug('Writing run report to run directory')
    _save_report(command_str, command_result,
                 start_time, entire_stdout, paths)
