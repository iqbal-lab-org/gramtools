import os
import logging

log = logging.getLogger('gramtools')


def _generate_project_paths(args):
    project_dir = args.gram_directory

    if hasattr(args, 'vcf') and args.vcf is not None:
        vcf_file_path = os.path.abspath(args.vcf)
    else:
        vcf_file_path = ''

    if hasattr(args, 'reference') and args.reference is not None:
        reference_file_path = os.path.abspath(args.reference)
    else:
        reference_file_path = ''

    def project_path(file_name):
        return os.path.join(project_dir, file_name)

    paths = {
        'project': project_dir,
        'vcf': vcf_file_path,
        'reference': reference_file_path,

        'prg': project_path('prg'),
        'encoded_prg': project_path('encoded_prg'),

        'variant_site_mask': project_path('variant_site_mask'),
        'allele_mask': project_path('allele_mask'),

        'fm_index': project_path('fm_index'),

        'perl_generated_vcf': project_path('perl_generated_vcf'),
        'perl_generated_fa': project_path('perl_generated_fa'),
    }
    return paths


def _move_perl_vcf(build_paths):
    original_fpath = build_paths['prg'] + '.vcf'
    target_fpath = build_paths['perl_generated_vcf']
    os.rename(original_fpath, target_fpath)


def _move_perl_reference(build_paths):
    original_fpath = build_paths['prg'] + '.fa'
    target_fpath = build_paths['perl_generated_fa']
    os.rename(original_fpath, target_fpath)


def perl_script_file_cleanup(build_paths):
    os.remove(build_paths['prg'] + '.mask_sites')
    os.remove(build_paths['prg'] + '.mask_alleles')
    _move_perl_vcf(build_paths)
    _move_perl_reference(build_paths)


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


def _quasimap_output_dirpath(outputs_base_path,
                             kmer_size,
                             start_time):
    template = '{start_time}_ksize{kmer_size}'
    directory_name = template.format(start_time=start_time, kmer_size=kmer_size)
    path = os.path.join(outputs_base_path, directory_name)
    return path


def generate_quasimap_paths(args, start_time):
    project_paths = _generate_project_paths(args)

    if hasattr(args, 'output_directory') and args.output_directory is not None:
        outputs_base_path = os.path.abspath(os.path.join(args.output_directory, os.pardir))
        run_dirpath = args.output_directory
    else:
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
