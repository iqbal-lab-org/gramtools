import os
import logging

log = logging.getLogger('gramtools')


def _generate_project_paths(args):
    project_dir = args.gram_directory
    kmer_prefix_diffs_file_name = \
        'kmer_prefix_diffs_' + str(args.kmer_size)
    kmer_index_file_name = \
        'kmer_index_' + str(args.kmer_size)

    if hasattr(args, 'vcf'):
        vcf_file_path = os.path.abspath(args.vcf)
    else:
        vcf_file_path = ''

    if hasattr(args, 'reference'):
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
        'sdsl_memory_log': project_path('sdsl_memory_log'),

        'kmers_dir': project_path('kmers'),
        'kmer_prefix_diffs': os.path.join(project_dir,
                                          'kmers',
                                          kmer_prefix_diffs_file_name),
        'kmer_index': os.path.join(project_dir,
                                   'kmers',
                                   kmer_index_file_name),

        'callgrind_out': project_path('callgrind.out'),

        'perl_generated_vcf': project_path('perl_generated_vcf'),
        'perl_generated_fa': project_path('perl_generated_fa'),
    }
    return paths


def _move_perl_sites_mask(build_paths):
    original_fpath = build_paths['prg'] + '.mask_sites'
    target_fpath = build_paths['variant_site_mask']
    os.rename(original_fpath, target_fpath)


def _move_perl_allele_mask(build_paths):
    original_fpath = build_paths['prg'] + '.mask_alleles'
    target_fpath = build_paths['allele_mask']
    os.rename(original_fpath, target_fpath)


def _move_perl_vcf(build_paths):
    original_fpath = build_paths['prg'] + '.vcf'
    target_fpath = build_paths['perl_generated_vcf']
    os.rename(original_fpath, target_fpath)


def _move_perl_reference(build_paths):
    original_fpath = build_paths['prg'] + '.fa'
    target_fpath = build_paths['perl_generated_fa']
    os.rename(original_fpath, target_fpath)


def perl_script_file_cleanup(build_paths):
    _move_perl_sites_mask(build_paths)
    _move_perl_allele_mask(build_paths)
    _move_perl_vcf(build_paths)
    _move_perl_reference(build_paths)


def generate_build_paths(args):
    project_paths = _generate_project_paths(args)
    return project_paths


def _generate_quasimap_run_dirpath(quasimap_outputs_dirpath,
                                   kmer_size,
                                   start_time):
    template = '{start_time}_ksize{kmer_size}'
    run_dirname = template.format(start_time=start_time,
                                  kmer_size=kmer_size)
    run_dirpath = os.path.join(quasimap_outputs_dirpath,
                               run_dirname)
    return run_dirpath


def generate_quasimap_paths(args, start_time):
    project_paths = _generate_project_paths(args)
    outputs_dirpath = os.path.join(project_paths['project'],
                                   'quasimap_outputs')
    run_dirpath = _generate_quasimap_run_dirpath(outputs_dirpath,
                                                 args.kmer_size,
                                                 start_time)

    def run_path(file_name):
        return os.path.join(run_dirpath, file_name)

    run_output_paths = {
        'quasimap_outputs_dirpath': outputs_dirpath,
        'quasimap_run_dirpath': run_dirpath,

        'run_report': run_path('report.json'),
        'callgrind_out': run_path('callgrind.out'),
        'sdsl_memory_log': run_path('sdsl_memory_log'),
        'allele_coverage': run_path('allele_coverage'),
        'reads_progress': run_path('reads_progress'),
    }

    paths = {
        'reads': args.reads
    }
    paths.update(project_paths)
    paths.update(run_output_paths)
    return paths


def check_project_file_structure(paths):
    log.debug('Checking project file structure')
    dir_names = [
        'project',
        'kmers_dir',
        'quasimap_outputs_dirpath',
        'quasimap_run_dirpath',
    ]
    for dir_name in dir_names:
        try:
            dir_path = paths[dir_name]
        except KeyError:
            continue
        if os.path.isdir(dir_path):
            continue
        log.debug('Creating missing directory:\n%s', dir_path)
        os.mkdir(dir_path)
