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


def get_run_dirpath(output_dirpath, species_dirpath, ksize):
    species_dir = os.path.basename(species_dirpath)
    time_str = str(time.time()).split('.')[0]

    template = '{time_str}_{species_dir}_ksize{ksize}'
    run_dir = template.format(time_str=time_str, species_dir=species_dir,
                              ksize=ksize)

    run_dirpath = os.path.join(output_dirpath, run_dir)
    return run_dirpath


def get_paths(args):
    # species_dirpath = get_species_dirpath(args.prg)
    species_dirpath = os.path.abspath(args.gram_files)
    output_dirpath = os.path.abspath(args.gram_files) + '_output'

    run_dirpath = get_run_dirpath(output_dirpath, species_dirpath, args.ksize)

    paths = {
        'species': species_dirpath,
        'prg': os.path.join(species_dirpath, 'prg'),
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

        'output': output_dirpath,
        'run': run_dirpath,
        'info': os.path.join(run_dirpath, 'info'),
        'fm_index_memory_log': os.path.join(
            run_dirpath, 'fm_index_memory_log'),
        'allele_coverage': os.path.join(run_dirpath, 'allele_coverage'),
        'reads': os.path.join(run_dirpath, 'reads'),
    }
    return paths


def setup_file_structure(paths):
    dirs = [
        paths['species'],
        paths['cache'],
        paths['kmer'],
        paths['kmer_suffix_array'],
        paths['output'],
        paths['run'],
    ]
    for dirpath in dirs:
        if os.path.isdir(dirpath):
            continue
        log.debug('Creating directory: %s', dirpath)
        os.mkdir(dirpath)


def execute_command(paths, args):
    command = [
        utils.gramtools_exec_fpath,
        '--prg', paths['prg'],
        '--csa', paths['fm_index'],
        '--ps', paths['sites_mask'],
        '--pa', paths['allele_mask'],
        '--co', paths['allele_coverage'],
        '--ro', paths['reads'],
        '--po', paths['int_encoded_prg'],
        '--log', paths['fm_index_memory_log'],
        '--kfile', paths['kmer_file'],
        '--input', paths['fast'],
        '--ksize', str(args.ksize),
    ]

    log.debug('Executing command:\n%s\n', ' '.join(command))

    current_working_directory = os.getcwd()
    process_handle = subprocess.Popen(command,
                                      cwd=current_working_directory,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
    utils.handle_process_result(process_handle)


def run(args):
    log.info('Start process: infer')

    paths = get_paths(args)
    setup_file_structure(paths)

    execute_command(paths, args)
    log.info('End process: infer')
