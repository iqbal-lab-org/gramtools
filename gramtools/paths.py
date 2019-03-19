## @file
# Sets up file names and directory structures for storing the gramtools inputs and outputs.
import os
import logging

log = logging.getLogger('gramtools')

## Defines the file names for all gramtools-related data files.
# These data files are committed to disk by the `build` process, and some/all used in later commands.
def _generate_project_paths(gram_dir):
    # if hasattr(args, 'reference') and args.reference is not None:
    #     reference_file_path = os.path.abspath(args.reference)
    # else:
    #     reference_file_path = ''

    def gram_path(file_name):
        return os.path.join(gram_dir, file_name)

    paths = {
        'gram_dir': gram_dir,
        # 'reference': reference_file_path,

        'prg': gram_path('prg'),
        'encoded_prg': gram_path('encoded_prg'),

        'variant_site_mask': gram_path('variant_site_mask'),
        'allele_mask': gram_path('allele_mask'),

        'fm_index': gram_path('fm_index'),

        'perl_generated': gram_path('perl_generated'),
        'perl_generated_vcf': gram_path('perl_generated.vcf'),
        'perl_generated_fa': gram_path('perl_generated.fa'),

        'build_report': gram_path('build_report.json'),

    }
    return paths


## Removes unnecessary accessory files produced by perl utility.
def perl_script_file_cleanup(build_paths):
    os.remove(build_paths['perl_generated'] + '.mask_sites')
    os.remove(build_paths['perl_generated'] + '.mask_alleles')
    os.rename(build_paths['perl_generated'],build_paths['prg'])

def generate_build_paths(args):
    paths = _generate_project_paths(args.gram_dir)

    if not os.path.exists(paths['gram_dir']):
        os.mkdir(paths['gram_dir'])
        log.debug('Creating gram directory:\n%s', paths['gram_dir'])
    return paths

def _quasimap_output_paths(quasimap_dir):
    path = lambda fname : os.path.join(quasimap_dir, fname)

    output_paths = {
        'allele_base_coverage': path('allele_base_coverage.json'),
        'grouped_allele_counts_coverage': path('grouped_allele_counts_coverage.json'),
        'allele_sum_coverage': path('allele_sum_coverage'),
        'read_stats' : path('read_stats.json'),
    }

    return output_paths

## Generates paths that will be shared between 'run' commands.
def _generate_run_paths(run_dir):

   def path_fact(base_dir):
       '''
        Factory building functions which build paths from a base_dir.
       '''
       return lambda fname: os.path.join(base_dir,fname)

   # Quasimap paths
   run_path = path_fact(run_dir)

   paths = {
       'run_dir' : run_dir,
       'build_dir' : run_path('build_dir'), # This will be a symlink of the actual build dir
       'reads_dir' : run_path('reads'), # This will contain symlinks to the actual reads
       'quasimap_dir' : run_path('quasimap_outputs'),
       'infer_dir' : run_path('infer_outputs'),
       'discover_dir' : run_path('discover_outputs'),
   }

   quasimap_path = path_fact(paths['quasimap_dir'])

   paths.update({

       'allele_base_coverage': quasimap_path('allele_base_coverage.json'),
       'grouped_allele_counts_coverage': quasimap_path('grouped_allele_counts_coverage.json'),
       'allele_sum_coverage': quasimap_path('allele_sum_coverage'),
       'read_stats' : quasimap_path('read_stats.json'),
   })

   # Infer paths
   infer_path = path_fact(paths['infer_dir'])

   paths.update({

       'inferred_fasta' : infer_path('inferred.fasta'),
       'inferred_vcf' : infer_path('inferred.vcf'),
       'inferred_ref_size' : infer_path('inferred_ref_size'),
   })

   # Discover paths
   discover_path = path_fact(paths['discover_dir'])
   paths.update({
       'cortex_vcf' : discover_path('cortex.vcf'),
       'rebased_vcf' : discover_path('rebased.vcf'),
   })

   return paths


## Make quasimap-related file and directory paths.
def generate_quasimap_paths(args):
    all_paths = _generate_project_paths(args)
    all_paths.update(_generate_run_paths(args.run_dir))

    if not os.path.exists(all_paths['run_dir']):
        os.mkdir(all_paths['run_dir'])
    else:
        log.warning("Directory {} already exists. Existing files and directories with conflicting names"
                    "(eg: 'reads', 'quasimap_outputs', 'build_dir') will be overriden.".format(all_paths['run_dir']))

    if not os.path.exists(all_paths['quasimap_dir']):
        os.mkdir(all_paths['quasimap_dir'])

    # Make a reference to the gram_dir made in build. This will avoid downstream commands
    # to require a --gram_dir argument.
    os.symlink(all_paths['gram_dir'], all_paths['build_dir'])

    if not os.path.exists(all_paths['reads_dir']):
        os.mkdir(all_paths['reads_dir'])

    # Make symlinks to the read files
    read_files = [read_file for arglist in args.reads for read_file in arglist], # Flattens out list of lists.
    for read_file in read_files:
        base = os.path.basename(read_file)
        os.symlink(read_file, os.path.join(all_paths['reads_dir'], base))

    all_paths.update({
        'report': os.path.join(all_paths['quasimap_dir'],'quasimap_report.json'),
    })

    return all_paths


def generate_infer_paths(args):

    all_paths = _generate_run_paths(args.run_dir)
    build_dir = all_paths['build_dir'] # Symlink to original --gram-dir, index against which quasimap was ran.
    all_paths.update(_generate_project_paths(build_dir))

    if not os.path.exists(all_paths['infer_dir']):
        os.mkdir(all_paths['infer_dir'])

    return all_paths


def generate_discover_paths(args):

    all_paths = _generate_run_paths(args.run_dir)
    build_dir = all_paths['build_dir'] # Symlink to original --gram-dir, index against which quasimap and infer ran.
    all_paths.update(_generate_project_paths(build_dir))

    if not os.path.exists(all_paths["discover_dir"]):
        os.mkdir(all_paths['discover_dir'])


    if not os.path.exists(all_paths['inferred_fasta']):
        log.error("Cannot find fasta formatted inferred personalised reference, at {}".format(_paths["inferred_fasta"]))

    if not os.path.exists(all_paths['inferred_vcf']):
        log.error("Cannot find vcf formatted inferred personalised reference, at {}".format(_paths["inferred_fasta"]))

    return all_paths



