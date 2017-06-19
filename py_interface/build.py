import os
import time
import logging
import argparse
import subprocess

from . import utils


log = logging.getLogger('gramtools')


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


def get_paths(args):
    species_dirpath = get_species_dirpath(args.vcf)

    paths = {
        'species': species_dirpath,
        'prg': os.path.join(species_dirpath, 'prg'),
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

        'perl_generated_fa': os.path.join(species_dirpath,
                                          'cache',
                                          'perl_generated_fa'),
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


def execute_command_generate_prg(paths, args):
    command = [
        'perl', utils.prg_build_exec_fpath,
        '--outfile', paths['prg'],
        '--vcf', paths['vcf'],
        '--ref', paths['fast'],
    ]

    log.debug('Executing command:\n%s\n', ' '.join(command))

    current_working_directory = os.getcwd()
    process_handle = subprocess.Popen(command,
                                      cwd=current_working_directory,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
    utils.handle_process_result(process_handle)
    log.debug('Finished executing command')


def execute_command_generate_kmers(paths, args):
    command = [
        'python2.7', utils.kmers_script_fpath,
        '-f', paths['perl_generated_fa'],
        '-k', str(args.ksize),
        '-n',
    ]

    log.debug('Executing command:\n%s\n', ' '.join(command))

    current_working_directory = os.getcwd()
    with open(paths['kmer_file'], 'wb') as kmers_fhandle:
        process_handle = subprocess.Popen(command,
                                          cwd=current_working_directory,
                                          stdout=kmers_fhandle,
                                          stderr=subprocess.PIPE)
    utils.handle_process_result(process_handle)
    log.debug('Finished executing command')


def file_cleanup_generate_prg(paths):
    original_fpath = paths['prg'] + '.mask_alleles'
    target_fpath = os.path.join(paths['species'], 'allele_mask')
    os.rename(original_fpath, target_fpath)

    original_fpath = paths['prg'] + '.mask_sites'
    target_fpath = os.path.join(paths['species'], 'sites_mask')
    os.rename(original_fpath, target_fpath)

    original_fpath = paths['prg']
    target_fpath = os.path.join(paths['species'], 'prg')
    os.rename(original_fpath, target_fpath)

    # TODO: should .prg.vcf be generated at all?
    original_fpath = paths['prg'] + '.vcf'
    target_fpath = os.path.join(paths['species'], 'cache',
                                'perl_generated_vcf')
    os.rename(original_fpath, target_fpath)

    original_fpath = paths['prg'] + '.fa'
    target_fpath = os.path.join(paths['species'], 'cache',
                                'perl_generated_fa')
    os.rename(original_fpath, target_fpath)


def run(args):
    log.info('Start process: build')

    paths = get_paths(args)
    setup_file_structure(paths)

    execute_command_generate_prg(paths, args)
    file_cleanup_generate_prg(paths)

    execute_command_generate_kmers(paths, args)

    log.info('End process: build')
