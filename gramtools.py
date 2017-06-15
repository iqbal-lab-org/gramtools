import os
import argparse
import logging


log = logging.getLogger('gramtools')
handler = logging.StreamHandler()
formatter = logging.Formatter(
    '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
log.addHandler(handler)
log.setLevel(logging.INFO)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("infer", help="",
                        action="store_true")
    parser.add_argument("prg", help="",
                        type=str)
    parser.add_argument("fast", help="",
                        type=str)
    parser.add_argument("ksize", help="",
                        type=int)
    parser.add_argument("output", help="",
                        type=int)
    args = parser.parse_args()
    return args


def check_path_exist(paths):
    missing_paths = set(path for path in paths
                        if not os.path.isfile(path))

    all_paths_present = not missing_paths
    if all_paths_present:
        return

    log.error("The following paths do not exist:")
    for path in missing_paths:
        log.error(path)
    exit(-1)


def get_species_dirpath(prg_fpath):
    prg_fname = os.path.basename(prg_fpath)
    if prg_fname.endswith('.prg'):
        species_dir = prg_fname[:-len('.prg')]
    else:
        species_dir = prg_fname

    current_working_directory = os.getcwd()
    species_dirpath = os.path.join(current_working_directory, species_dir)
    return species_dirpath


def get_prg_fpath(species_dirpath, prg_arg_fpath):
    prg_fname = os.path.basename(prg_arg_fpath)
    prg_fpath = os.path.join(species_dirpath, prg_fname)
    return prg_fpath


def get_run_dirpath(species_dirpath, fast_fpath, ksize):
    species_dir = species_dirpath.split('/')[-1]

    template = 'run_{species_dir}_{fast}_{ksize}'
    run_dir = template.format(species_dir=species_dirpath,
                              fast=fast_fpath, ksize=ksize)

    run_dirpath = os.path.join(output_dirpath, run_dir)
    return run_dirpath


def get_paths(args):
    species_dirpath = get_species_dirpath(args.prg)
    run_dirpath = get_run_dirpath(species_dirpath, arg.fast, args.ksize)

    paths = {
        'species': species_dirpath,
        'prg': os.path.join(species_dirpath, 'prg'),
        'sites_mask': os.path.join(species_dirpath, 'sites_mask'),
        'allele_mask': os.path.join(species_dirpath, 'allele_mask'),
        'fast': arg.fast,

        'kmer': os.path.join(species_dirpath, 'kmer'),
        'cache': os.path.join(species_dirpath, 'cache'),
        'int_encoded_prg': os.path.join(species_dirpath, 'cache', 'int_encoded_prg'),
        'fm_index': os.path.join(species_dirpath, 'cache', 'fm_index'),
        'kmer_suffix_array': os.path.join(species_dirpath, 'cache',
                                          'kmer_suffix_array'),

        'output': args.output,
        'run': run_dirpath,
        'info': os.path.join(run_dirpath, 'info'),
        'fm_index_memory_log': os.path.join(run_dirpath, 'fm_index_memory_log'),
        'allele_coverage': os.path.join(run_dirpath, 'allele_coverage'),
        'processed_reads': os.path.join(run_dirpath, 'processed_reads'),
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
    new_dirs = [dir for dir in dirs if not os.path.isdir(dir)]
    map(os.mkdir, new_dirs)


def populate_file_structure(paths, args):
    if not os.path.isfile(paths['prg']):
        os.symlink(arg.prg, paths['prg'])


def generate_sites_mask():
    pass


def infer(args):
    check_path_exist([arg.prg, arg.fast], os.path.isfile)
    paths = get_paths(args)

    setup_file_structure(paths)
    populate_file_structure(paths, args)

    generate_sites_mask()


if __name__ == '__main__'
    args = parse_args()

    if args.infer:
        infer(args)