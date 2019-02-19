## @file
# Sets up file names and directory structures for storing the gramtools inputs and outputs.
import os
import uuid
import logging
import tempfile

log = logging.getLogger('gramtools')

## Defines the file names for all gramtools-related data files.
# These data files are committed to disk by the `build` process, and loaded for the `quasimap` process.
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
    }
    return paths


## Removes unnecessary accessory files produced by perl utility.
def perl_script_file_cleanup(build_paths):
    os.remove(build_paths['prg'] + '.mask_sites')
    os.remove(build_paths['prg'] + '.mask_alleles')

def generate_build_paths(args):
    project_paths = _generate_project_paths(args)
    project_paths['build_report'] = os.path.join(project_paths['project'], 'build_report.json')
    return project_paths


def generate_quasimap_run_paths(args):
    def path(fname):
        return os.path.join(args.quasimap_directory, fname)

    return {
        'allele_base_coverage': path('allele_base_coverage.json'),
        'grouped_allele_counts_coverage': path('grouped_allele_counts_coverage.json'),
        'allele_sum_coverage': path('allele_sum_coverage'),
        'report': path('report.json'),
    }


def generate_infer_paths(args):
    project_paths = _generate_project_paths(args)
    quasimap_paths = generate_quasimap_run_paths(args)
    paths = {**project_paths, **quasimap_paths}
    return paths


def generate_discover_paths(args):
    def add_uuid(template: str):
        return template.format(uuid=str(uuid.uuid4()))

    project_paths = _generate_project_paths(args)
    tmp_base_directory = tempfile.gettempdir()
    tmp_directory_name = add_uuid('gramtools_discover_{uuid}')
    tmp_directory = os.path.join(tmp_base_directory, tmp_directory_name)

    base_reference_fasta_file_name = add_uuid('base_reference_{uuid}.fasta')
    inferred_fasta_file_name = add_uuid('infer_{uuid}.fasta')
    cortex_vcf_file_name = add_uuid('cortex_{uuid}.vcf')
    rebase_vcf_file_name = add_uuid('rebase_{uuid}.vcf')

    _paths = {
        **project_paths,
        'tmp_directory': tmp_directory,
        'base_reference': os.path.join(tmp_directory, base_reference_fasta_file_name),
        'inferred_reference': os.path.join(tmp_directory, inferred_fasta_file_name),
        'cortex_vcf': os.path.join(tmp_directory, cortex_vcf_file_name),
        'rebase_vcf': os.path.join(tmp_directory, rebase_vcf_file_name),
    }
    return _paths


def _quasimap_output_dirpath(outputs_base_path,
                             kmer_size,
                             start_time):
    template = '{start_time}_ksize{kmer_size}'
    directory_name = template.format(start_time=start_time, kmer_size=kmer_size)
    path = os.path.join(outputs_base_path, directory_name)
    return path

## Make quasimap-related file and directory paths.
def generate_quasimap_paths(args, start_time):
    project_paths = _generate_project_paths(args)

    if hasattr(args, 'output_directory') and args.output_directory is not None:
        outputs_base_path = os.path.abspath(os.path.join(args.output_directory, os.pardir))
        run_dirpath = args.output_directory
    else: # Name a default output directory inside the gram project directory.
        outputs_base_path = os.path.join(project_paths['project'], 'quasimap_outputs')
        run_dirpath = _quasimap_output_dirpath(outputs_base_path,
                                               args.kmer_size,
                                               start_time)

    paths = {
        'reads': args.reads,

        'quasimap_outputs_dirpath': outputs_base_path,
        'quasimap_run_dirpath': run_dirpath,
        'run_report': os.path.join(run_dirpath, 'report.json'),
    }
    paths.update(project_paths)
    return paths

## Makes project directories if they are required and do not exist.
def check_project_file_structure(paths):
    log.debug('Checking project file structure')
    directory_keys = [
        'project',
        'quasimap_outputs_dirpath',
        'quasimap_run_dirpath',
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
