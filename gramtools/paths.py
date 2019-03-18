## @file
# Sets up file names and directory structures for storing the gramtools inputs and outputs.
import os
import logging

log = logging.getLogger('gramtools')

## Defines the file names for all gramtools-related data files.
# These data files are committed to disk by the `build` process, and some/all used in later commands.
def _generate_project_paths(args):
    project_dir = args.gram_directory

    if hasattr(args, 'reference') and args.reference is not None:
        reference_file_path = os.path.abspath(args.reference)
    else:
        reference_file_path = ''

    def project_path(file_name):
        return os.path.join(project_dir, file_name)

    paths = {
        'project': project_dir,
        'reference': reference_file_path,

        'prg': project_path('prg'),
        'encoded_prg': project_path('encoded_prg'),

        'variant_site_mask': project_path('variant_site_mask'),
        'allele_mask': project_path('allele_mask'),

        'fm_index': project_path('fm_index'),

        'perl_generated': project_path('perl_generated'),
        'perl_generated_vcf': project_path('perl_generated.vcf'),
        'perl_generated_fa': project_path('perl_generated.fa'),

        'build_report': project_path('build_report.json')
    }
    return paths


## Removes unnecessary accessory files produced by perl utility.
def perl_script_file_cleanup(build_paths):
    os.remove(build_paths['perl_generated'] + '.mask_sites')
    os.remove(build_paths['perl_generated'] + '.mask_alleles')
    os.rename(build_paths['perl_generated'],build_paths['prg'])

def generate_build_paths(args):
    project_paths = _generate_project_paths(args)
    return project_paths


def generate_quasimap_run_paths(args):
    def path(fname):
        return os.path.join(args.quasimap_directory, fname)

    return {
        'allele_base_coverage': path('allele_base_coverage.json'),
        'grouped_allele_counts_coverage': path('grouped_allele_counts_coverage.json'),
        'allele_sum_coverage': path('allele_sum_coverage'),
        'report': path('report.json'),
        'read_stats' : path('read_stats.json'),
    }


def generate_infer_paths(args):
    project_paths = _generate_project_paths(args)
    quasimap_paths = generate_quasimap_run_paths(args)
    paths = {**project_paths, **quasimap_paths}

    if args.infer_dir is not None:
        paths["infer_dir"] = args.infer_dir
    else:
        paths["infer_dir"] = os.path.join(paths["project"],"infer_outputs")

    if not os.path.exists(paths["infer_dir"]):
        os.mkdir(paths["infer_dir"])

    paths["inferred_fasta"] = os.path.join(paths["infer_dir"],"inferred.fasta")
    paths["inferred_vcf"] = os.path.join(paths["infer_dir"],"inferred.vcf")
    paths["inferred_ref_size"] = os.path.join(paths["infer_dir"], "inferred_ref_size")

    return paths


def generate_discover_paths(args):

    project_paths = _generate_project_paths(args)

    _paths = {
        **project_paths,
        "infer_dir" : os.path.abspath(args.infer_dir),
    }

    # Discovery directory for permanent output
    if args.discover_dir is not None:
        _paths["discover_dir"] = args.discover_dir
    else:
        _paths["discover_dir"] = os.path.join(os.path.abspath(_paths["project"]), "discover_outputs")

    if not os.path.exists(_paths["discover_dir"]):
        os.mkdir(_paths["discover_dir"])

    # Outputs
    _paths['cortex_vcf'] = os.path.join(_paths["discover_dir"], "cortex.vcf")
    _paths['rebased_vcf'] = os.path.join(_paths["discover_dir"], "rebased.vcf")


    # Check we have the required files from `infer`.
    _paths["inferred_fasta"] = os.path.join(_paths["infer_dir"],"inferred.fasta")
    _paths["inferred_vcf"] = os.path.join(_paths["infer_dir"],"inferred.vcf")
    _paths["inferred_ref_size"] = os.path.join(_paths["infer_dir"], "inferred_ref_size")

    if not os.path.exists(_paths["inferred_fasta"]):
        log.error("Cannot find fasta formatted inferred personalised reference, at {}".format(_paths["inferred_fasta"]))

    if not os.path.exists(_paths["inferred_vcf"]):
        log.error("Cannot find vcf formatted inferred personalised reference, at {}".format(_paths["inferred_fasta"]))


    return _paths


## Make quasimap-related file and directory paths.
def generate_quasimap_paths(args, start_time):
    project_paths = _generate_project_paths(args)

    quasimap_dir = args.quasimap_dir
    if quasimap_dir is None: # Name a default output directory inside the gram project directory.
        quasimap_dir = os.path.join(project_paths['project'], 'quasimap_outputs')

    paths = {
        'reads': [read_file for arglist in args.reads for read_file in arglist], # Flattens out list of lists

        'quasimap_dir': quasimap_dir,
        'run_report': os.path.join(quasimap_dir, 'quasimap_report.json'),
    }
    paths.update(project_paths)
    return paths

## Makes project directories if they are required and do not exist.
def check_project_file_structure(paths):
    log.debug('Checking project file structure')
    directory_keys = [
        'project',
        'quasimap_dir',
    ]
    for directory_key in directory_keys:
        try:
            directory_path = paths[directory_key]
        except KeyError:
            continue
        if os.path.isdir(directory_path):
            continue
        log.debug('Creating missing directory:\n%s', directory_path)
        os.mkdir(directory_path)
