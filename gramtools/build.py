import os
import copy
import time
import logging
import subprocess

from . import common
from . import kmers

log = logging.getLogger('gramtools')


def parse_args(common_parser, subparsers):
    parser = subparsers.add_parser('build',
                                   parents=[common_parser])

    parser.add_argument('--vcf', help='',
                        type=str,
                        required=True)
    parser.add_argument('--reference', help='',
                        type=str,
                        required=True)
    parser.add_argument('--kmer-size', help='',
                        type=int,
                        required=True)
    parser.add_argument('--kmer-region-size',
                        dest='kmer_region_size',
                        help='',
                        type=int,
                        required=True)

    parser.add_argument('--nonvariant-kmers', help='',
                        default=False,
                        action='store_true')


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


def _get_paths(args):
    project_dirpath = _get_project_dirpath(args.vcf)

    paths = {
        'project': project_dirpath,
        'prg': os.path.join(project_dirpath, 'prg'),
        'vcf': os.path.abspath(args.vcf),
        'sites_mask': os.path.join(project_dirpath, 'sites_mask'),
        'allele_mask': os.path.join(project_dirpath, 'allele_mask'),
        'reference': args.reference,

        'kmer': os.path.join(project_dirpath, 'kmer'),
        'kmer_file': os.path.join(project_dirpath, 'kmer',
                                  'ksize_' + str(args.kmer_size)),
        'cache': os.path.join(project_dirpath, 'cache'),
        'int_encoded_prg': os.path.join(
            project_dirpath, 'cache', 'int_encoded_prg'),
        'fm_index': os.path.join(project_dirpath, 'cache', 'fm_index'),
        'kmer_suffix_array': os.path.join(project_dirpath, 'cache',
                                          'kmer_suffix_array'),

        'perl_generated_fa': os.path.join(project_dirpath,
                                          'cache',
                                          'perl_generated_fa'),
    }
    return paths


def _setup_file_structure(paths):
    dirs = [
        paths['project'],
        paths['cache'],
        paths['kmer'],
        paths['kmer_suffix_array'],
    ]
    for dirpath in dirs:
        if os.path.isdir(dirpath):
            continue
        log.debug('Creating directory: %s', dirpath)
        os.mkdir(dirpath)


def _execute_command_generate_prg(paths, _):
    command = [
        'perl', common.prg_build_exec_fpath,
        '--outfile', paths['prg'],
        '--vcf', paths['vcf'],
        '--ref', paths['reference'],
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
    log.debug('Finished executing command: %.3f seconds', timer_end - timer_start)


def _file_cleanup_generate_prg(paths):
    original_fpath = paths['prg'] + '.mask_alleles'
    target_fpath = os.path.join(paths['project'], 'allele_mask')
    os.rename(original_fpath, target_fpath)

    original_fpath = paths['prg'] + '.mask_sites'
    target_fpath = os.path.join(paths['project'], 'sites_mask')
    os.rename(original_fpath, target_fpath)

    original_fpath = paths['prg']
    target_fpath = os.path.join(paths['project'], 'prg')
    os.rename(original_fpath, target_fpath)

    # TODO: should .prg.vcf be generated at all?
    original_fpath = paths['prg'] + '.vcf'
    target_fpath = os.path.join(paths['project'], 'cache',
                                'perl_generated_vcf')
    os.rename(original_fpath, target_fpath)

    original_fpath = paths['prg'] + '.fa'
    target_fpath = os.path.join(paths['project'], 'cache',
                                'perl_generated_fa')
    os.rename(original_fpath, target_fpath)


def _generate_kmers(paths, args):
    log.debug('Generating kmers from PRG')
    timer_start = time.time()

    args = copy.copy(args)
    args.reference = paths['perl_generated_fa']
    args.output_fpath = paths['kmer_file']
    args.sites_mask_fpath = paths['sites_mask']
    args.allele_mask_fpath = paths['allele_mask']

    kmers.run(args)

    timer_end = time.time()
    log.debug('Finished executing command: %.3f seconds', timer_end - timer_start)


def run(args):
    log.info('Start process: build')

    paths = _get_paths(args)
    _setup_file_structure(paths)

    _execute_command_generate_prg(paths, args)
    _file_cleanup_generate_prg(paths)

    _generate_kmers(paths, args)

    log.info('End process: build')
