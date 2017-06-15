import os
import subprocess

import sacred
from sacred.utils import apply_backspaces_and_linefeeds
from sacred.observers.file_storage import FileStorageObserver


GRAMTOOLS_INSTALL_PATH = './gramtools'
GENERATE_PRG_SCRIPT_PATH = os.path.join(GRAMTOOLS_INSTALL_PATH,
                                        'utils/vcf_to_linear_prg.pl')
GENERATE_KMERS_SCRIPT_PATH = os.path.join(GRAMTOOLS_INSTALL_PATH,
                                          'utils/variantKmers.py')
MAP_READS_PATH = os.path.join(GRAMTOOLS_INSTALL_PATH, 'bin',
                              'gramtools')


experiment = sacred.Experiment()
experiment.captured_out_filter = apply_backspaces_and_linefeeds

file_observer = FileStorageObserver.create('gramtools_runs')
file_observer.save_sources = lambda x: None
experiment.observers.append(file_observer)


def generate_paths(vcf_path, fasta_path):
    """Generate and return all file paths associated with experiment."""
    vcf_path = os.path.abspath(vcf_path)
    file_observer.run_entry['artifacts'].append(vcf_path)

    fasta_path = os.path.abspath(fasta_path)
    file_observer.run_entry['artifacts'].append(fasta_path)

    run_path = os.path.abspath(file_observer.dir)
    data_path = os.path.join(run_path, 'data')

    def make_path(fname): return os.path.join(data_path, fname)

    paths = {
        'run_path': run_path,
        'data_path': data_path,
        'vcf': vcf_path,
        'fasta': fasta_path,

        'modified_fasta': make_path('modified_fasta'),
        'prg': make_path('prg'),
        'kmers': make_path('kmers'),
        'mask_sites': make_path('mask_sites'),
        'mask_alleles': make_path('mask_alleles'),
        'coverage': make_path('coverage'),
        'csa': make_path('csa'),
        'csa_memory_log': make_path('csa_memory_log'),
        'reads': make_path('reads'),
        'prg_integer_encoding': make_path('prg_integer_encoding'),
    }
    return paths


def handle_process_result(process_handle, _log):
    """Report process results to logging."""
    stdout, error_message = process_handle.communicate()
    error_code = process_handle.returncode

    if error_message:
        _log.info('Process error message:\n%s',
                  error_message.decode("utf-8"))
        _log.info('Process error code: %s', error_code)

    if stdout:
        truncate_stdout = stdout.decode("utf-8")[:100]
        _log.info('Truncated stdout:\n%s', truncate_stdout)
    else:
        _log.info('stdout: none or piped out')

    if error_code != 0:
        raise RuntimeError('Error code != 0')


@experiment.capture
def generate_prg(paths, _log):
    """Generate a Population Reference Graph (PRG) file."""
    process_name = os.path.basename(GENERATE_PRG_SCRIPT_PATH)
    _log.info('Start process: %s', process_name)

    command = ['perl', GENERATE_PRG_SCRIPT_PATH,
               '--outfile', paths['prg'],
               '--vcf', paths['vcf'],
               '--ref', paths['fasta']]
    process_handle = subprocess.Popen(command, cwd=paths['run_path'],
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
    handle_process_result(process_handle, _log)

    rename_mapping = [
        ('prg.fa', 'modified_fasta'),
        ('prg.mask_alleles', 'mask_alleles'),
        ('prg.mask_sites', 'mask_sites'),
        ('prg.vcf', 'vcf'),
    ]

    for old_fname, new_fname in rename_mapping:
        old_path = os.path.join(paths['data_path'], old_fname)
        new_path = os.path.join(paths['data_path'], new_fname)
        os.rename(old_path, new_path)

    _log.info('End process: %s', process_name)


@experiment.capture
def generate_kmers(kmer_size, paths, _log):
    """Generate K-mers file."""
    process_name = os.path.basename(GENERATE_KMERS_SCRIPT_PATH)
    _log.info('Start process: %s', process_name)

    command = ['python2.7', GENERATE_KMERS_SCRIPT_PATH,
               '-f', paths['modified_fasta'],
               '-k', str(kmer_size),
               '-n']
    with open(paths['kmers'], 'wb') as kmers_fhandle:
        process_handle = subprocess.Popen(command, cwd=paths['run_path'],
                                          stdout=kmers_fhandle,
                                          stderr=subprocess.PIPE)
    handle_process_result(process_handle, _log)
    _log.info('End process: %s', process_name)


@experiment.capture
def map_reads_to_prg(kmer_size, paths, _log):
    """Map reads to PRG."""
    process_name = os.path.basename(MAP_READS_PATH)
    _log.info('Start process: %s', process_name)

    command = [MAP_READS_PATH,
               '--prg', paths['prg'],
               '--csa', paths['csa'],
               '--ps', paths['mask_sites'],
               '--pa', paths['mask_alleles'],
               '--co', paths['coverage'],
               '--ro', paths['reads'],
               '--po', paths['prg_integer_encoding'],
               '--log', paths['csa_memory_log'],
               '--kfile', paths['kmers'],
               '--input', paths['fasta'],
               '--ksize', str(kmer_size)]

    process_handle = subprocess.Popen(command, cwd=paths['run_path'],
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
    handle_process_result(process_handle, _log)
    _log.info('End process: %s', process_name)


@experiment.config
def config():
    """Configuration paramaters for experiment manager."""
    vcf_path = ''
    fasta_path = ''
    kmer_size = 9


@experiment.automain
def experiment_manager(_log, vcf_path='', fasta_path='', kmer_size=9):
    """Run experiments by executing tool chain."""
    paths = generate_paths(vcf_path, fasta_path)
    os.mkdir(os.path.join(paths['run_path'], 'data'))

    try:
        generate_prg(paths, _log)
        generate_kmers(kmer_size, paths, _log)
        map_reads_to_prg(kmer_size, paths, _log)
    except RuntimeError as e:
        _log.error(e)
        return e
