import os
import time
import argparse
import logging
import subprocess

import utils
from utils import log, prg_build_exec_fpath, kmers_script_fpath


def get_species_dirpath(prg_fpath):
    prg_fname = os.path.basename(prg_fpath)
    if prg_fname.endswith('.prg'):
        species_dir = prg_fname[:-len('.prg')]
    else:
        species_dir = prg_fname
    species_dir = species_dir.split('.')[0]

    current_working_directory = os.getcwd()
    species_dirpath = os.path.join(current_working_directory, species_dir)
    return species_dirpath


def get_prg_fpath(species_dirpath, prg_arg_fpath):
    prg_fname = os.path.basename(prg_arg_fpath)
    prg_fpath = os.path.join(species_dirpath, prg_fname)
    return prg_fpath


def get_paths(args):
    species_dirpath = get_species_dirpath(args.prg)

    paths = {
        'species': species_dirpath,
        'prg': os.path.abspath(args.prg),
        'vcf': os.path.abspath(args.vcf),
        'sites_mask': os.path.join(species_dirpath, 'sites_mask'),
        'allele_mask': os.path.join(species_dirpath, 'allele_mask'),
        'fast': args.fast,

        'kmer': os.path.join(species_dirpath, 'kmer'),
        'kmer_file': os.path.join(species_dirpath, 'kmer',
                                  'ksize_' + str(args.ksize)),
        'cache': os.path.join(species_dirpath, 'cache'),
        'int_encoded_prg': os.path.join(
            species_dirpath, 'cache', 'int_encoded_prg'),
        'fm_index': os.path.join(species_dirpath, 'cache', 'fm_index'),
        'kmer_suffix_array': os.path.join(species_dirpath, 'cache',
                                          'kmer_suffix_array'),

        'output': args.output,
    }
    return paths


def setup_file_structure(paths):
    dirs = [
        paths['species'],
        paths['cache'],
        paths['kmer'],
        paths['kmer_suffix_array'],
    ]
    for dirpath in dirs:
        if os.path.isdir(dirpath):
            continue
        log.debug('Creating directory: %s', dirpath)
        os.mkdir(dirpath)


def execute_command(paths, args):
    command = [
        'perl', prg_build_exec_fpath,
        '--outfile', paths['prg'],
        '--vcf', paths['vcf'],
        '--ref', paths['fast'],
    ]

    log.debug('Executing command:\n%s', '\n'.join(command))

    current_working_directory = os.getcwd()
    process_handle = subprocess.Popen(command,
                                      cwd=current_working_directory,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
    utils.handle_process_result(process_handle)

    command = [
        'python2.7', kmers_script_fpath,
        '-f', paths['modified_fasta'],
        '-k', str(kmer_size),
        '-n',
    ]
    with open(paths['kmer'], 'wb') as kmers_fhandle:
        process_handle = subprocess.Popen(command, cwd=paths['run_path'],
                                          stdout=kmers_fhandle,
                                          stderr=subprocess.PIPE)
    utils.handle_process_result(process_handle)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--build", help="",
                        action="store_true")

    parser.add_argument("--prg", help="",
                        type=str)
    parser.add_argument("--fast", help="",
                        type=str)
    parser.add_argument("--vcf", help="",
                        type=str)

    args = parser.parse_args()
    return args


def file_cleanup(paths):
    original_fpath = paths['prg'] + '.mask_alleles'
    target_fpath = os.path.join(paths['species'], 'allele_mask')
    os.rename(original_fpath, target_fpath)

    original_fpath = paths['prg'] + '.mask_sites'
    target_fpath = os.path.join(paths['species'], 'sites_mask')
    os.rename(original_fpath, target_fpath)


def run(args):
    log.info('Start process: build')

    utils.check_path_exist([args.vcf, args.fast])
    paths = get_paths(args)
    setup_file_structure(paths)

    execute_command(paths, args)
    file_cleanup(paths)

    log.info('End process: build')
