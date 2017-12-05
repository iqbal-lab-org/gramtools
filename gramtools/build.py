import os
import time
import logging
import subprocess

from . import common
from . import paths

log = logging.getLogger('gramtools')


def parse_args(common_parser, subparsers):
    parser = subparsers.add_parser('build',
                                   parents=[common_parser])
    parser.add_argument('--gram-directory',
                        help='',
                        type=str,
                        required=True)
    parser.add_argument('--vcf',
                        help='',
                        type=str,
                        required=True)
    parser.add_argument('--reference',
                        help='',
                        type=str,
                        required=True)

    parser.add_argument('--kmer-size',
                        help='',
                        type=int,
                        default=15,
                        required=False)
    parser.add_argument('--kmer-region-size',
                        help='',
                        type=int,
                        default=150,
                        required=False)
    parser.add_argument('--max-read-length',
                        help='',
                        type=int,
                        default=150,
                        required=False)

    parser.add_argument('--all-kmers',
                        help='',
                        action='store_true',
                        required=False)


def _execute_command_generate_prg(build_paths, _):
    command = [
        'perl', common.prg_build_exec_fpath,
        '--outfile', build_paths['prg'],
        '--vcf', build_paths['vcf'],
        '--ref', build_paths['reference'],
    ]

    log.debug('Executing command:\n\n%s\n', ' '.join(command))
    timer_start = time.time()

    current_working_directory = os.getcwd()
    process_handle = subprocess.Popen(command,
                                      cwd=current_working_directory,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
    common.handle_process_result(process_handle)
    timer_end = time.time()
    log.debug('Finished executing command: %.3f seconds',
              timer_end - timer_start)


def _execute_gramtools_cpp_build(build_paths, args):
    command = [
        common.gramtools_exec_fpath,
        'build',
        '--gram', build_paths['project'],
        '--prg', build_paths['prg'],
        '--encoded-prg', build_paths['encoded_prg'],
        '--fm-index', build_paths['fm_index'],
        '--variant-site-mask', build_paths['variant_site_mask'],
        '--allele-mask', build_paths['allele_mask'],
        '--memory-log', build_paths['sdsl_memory_log'],
        '--kmer-index', build_paths['kmer_index'],
        '--kmer-size', str(args.kmer_size),
        '--max-read-size', str(args.max_read_length),
    ]

    if args.debug:
        command += ['--debug']

    callgrind_command = [
        'valgrind',
        '--tool=callgrind',
        '--callgrind-out-file=' + build_paths['callgrind_out'],
    ]

    callgrind_command = callgrind_command if args.profile else []
    command = callgrind_command + command
    command_str = ' '.join(command)

    log.debug('Executing command:\n\n%s\n', command_str)

    current_working_directory = os.getcwd()
    log.debug('Using current working directory:\n%s',
              current_working_directory)

    process_handle = subprocess.Popen(command_str,
                                      cwd=current_working_directory,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE,
                                      shell=True)

    process_result = common.handle_process_result(process_handle)
    command_result, entire_stdout = process_result
    return command_str, command_result, entire_stdout


def run(args):
    log.info('Start process: build')
    if hasattr(args, 'max_read_length'):
        args.kmer_region_size = args.max_read_length

    build_paths = paths.generate_build_paths(args)
    paths.check_project_file_structure(build_paths)

    _execute_command_generate_prg(build_paths, args)
    paths.perl_script_file_cleanup(build_paths)

    _execute_gramtools_cpp_build(build_paths, args)

    log.info('End process: build')
