## @file
#  Build/load a population reference genome and set it up for quasimapping.
# Either a vcf/reference is passed and a prg generated from it, or an existing prg is passed.
# Once the prg is stored the back-end `build` routine is called, producing the encoded prg, its fm-index, and other supporting data structures.
import os
import time
import json
import shutil
import logging
import subprocess
import collections

import cluster_vcf_records

from .. import version
from .. import common
from .. import paths
from ..utils import vcf_to_prg_string

log = logging.getLogger('gramtools')


def parse_args(common_parser, subparsers):
    parser = subparsers.add_parser('build',
                                   parents=[common_parser])
    parser.add_argument('--gram-dir', '--gram-directory',
                        help='',
                        dest='gram_dir',
                        type=str,
                        required=True)

    parser.add_argument('--vcf',
                        help='File containing variant information to capture in the prg.',
                        nargs="+",
                        action="append",
                        type=str)
    parser.add_argument('--reference',
                        help='Reference the vcf file refers to, used to build non-variant parts of the prg.',
                        type=str,
                        required=False)
    parser.add_argument('--prg',
                        help='A prg can be passed in directly instead of a vcf/reference combination.',
                        type=str,
                        required=False)

    parser.add_argument('--kmer-size',
                        help='Kmer size for indexing the prg. Defaults to 5.',
                        type=int,
                        default=5,
                        required=False)


    # The current default behaviour is to extract only relevant kmers from prg.
    parser.add_argument('--all-kmers',
                        help='Whether or not all kmers of given size should be indexed.\n'
                             'When this flag is not used, only kmers overlapping variant sites in prg will be indexed.',
                        action='store_true',
                        required=False)

    parser.add_argument('--max-read-length',
                        help='Used to determine which kmers overlap variant sites. Only needed if --all-kmers flag is off.',
                        type=int,
                        default=150,
                        required=False)

    parser.add_argument('--max-threads',
                        help='',
                        type=int,
                        default=1,
                        required=False)


def with_report(f):
    """
    Decorator to add logging and reporting to build procedures
    To signal that something went wrong in decorated function call, that function needs to raise an Exception.
    """
    def reportify(report, action, *args):
        if report.get('success') is False:
            report[action] = {
                'success': False
            }
            return report

        success, error, command, stdout = True, None, None, None
        timer_start = time.time()

        try:
            original_result = f(report, action, *args)
        except Exception as e:
            success = False
            error = str(e)
            original_result = None
        timer_end = time.time()

        log.debug(f'Ran {action} in: {timer_end - timer_start} seconds')

        report['success'] = success
        action_report = collections.OrderedDict([
        ('success', success),
        ('error_message', error),
            ('Run time', int(timer_end) - int(timer_start))
        ])

        # The condition below allows the called function to modify the report as well.
        if action not in report:
            report[action] = action_report
        else:
            report[action].update(action_report)

        return original_result
    return reportify


def run(args):
    log.info('Start process: build')
    start_time = str(time.time()).split('.')[0]

    command_paths = paths.generate_build_paths(args)
    report = collections.OrderedDict()

    # Update the vcf path to a combined vcf from all those provided.
    # We also do this if only a single one is provided, to deal with overlapping records.
    command_paths['vcf'] = _cluster_vcf_records(report, 'vcf_record_clustering', command_paths)

    if hasattr(args, 'prg') and args.prg is not None:
        _skip_prg_construction(report, 'copy_existing_PRG_string', command_paths, args)
    else:
        _execute_command_generate_prg(report, 'vcf_to_PRG_string_conversion', command_paths)

    _execute_gramtools_cpp_build(report, 'gramtools_build', command_paths, args)

    log.debug('Computing sha256 hash of project paths')
    command_hash_paths = common.hash_command_paths(command_paths)

    log.debug('Saving command report:\n%s', command_paths['build_report'])
    _report = collections.OrderedDict([
        ('kmer_size', args.kmer_size),
        ('max_read_length', args.max_read_length)
    ])
    report.update(_report)

    report_file_path = command_paths['build_report']
    _save_report(start_time,
                 report,
                 command_paths,
                 command_hash_paths,
                 report_file_path)

    if report["success"] is False:
        log.error(f'Unsuccessful build. Process reported to {command_paths["build_report"]}')
        exit(1)
    else:
        log.info(f'Success! Build process report in {command_paths["build_report"]}')


def _count_vcf_record_lines(vcf_file_path):
    num_recs = 0
    with open(vcf_file_path) as f_in:
        for line in f_in:
            if line[0] != "#":
                num_recs += 1
    return num_recs


## Combines records in one or more vcf files using external python utility.
# Records where the REFs overlap are merged together and all possible haplotypes enumerated.
# New path to 'vcf' is path to the combined vcf
@with_report
def _cluster_vcf_records(report, action, command_paths):
    vcf_files = command_paths['vcf']
    final_vcf_name = command_paths['built_vcf']
    log.info(f"Running {action} on {str(vcf_files)}.")

    cluster = cluster_vcf_records.vcf_clusterer.VcfClusterer(vcf_files,
                                                             command_paths['original_reference'],
                                                             final_vcf_name, max_alleles_per_cluster=5000)
    cluster.run()

    return final_vcf_name


## Calls utility that converts a vcf and fasta reference into a linear prg.
@with_report
def _execute_command_generate_prg(report, action, build_paths):
    vcf_in = build_paths['vcf']
    log.info(f"Running {action} on {vcf_in}")

    converter = vcf_to_prg_string.Vcf_to_prg(vcf_in, build_paths['original_reference'],
                                             build_paths['prg_string'], mode="normal")
    # print(converter.prg_vector)
    converter._write_bytes()

    ## The converter does not produce a vcf
    # Thus the input vcf needs to be 'clean': we need each vcf record to be converted
    # to a variant site in the prg so that later on the genotyping process, which uses vcf, is possible.
    num_recs_in_vcf = _count_vcf_record_lines(vcf_in)
    assert num_recs_in_vcf == converter.num_sites, log.error(f"Mismatch between number of vcf records in {vcf_in}"
                                                             f"({num_recs_in_vcf} and number of variant sites in"
                                                             f"PRG string ({converter.num_sites}.\n"
                                                             f"Possible source of error: vcf record clustering does not"
                                                             f"produce non-overlapping records, or conversion utility"
                                                             f" {vcf_to_prg_string.__file__} is broken.")

## Checks prg file exists and copies it to gram directory.
@with_report
def _skip_prg_construction(report, action, build_paths, args):
    if not os.path.isfile(args.prg):
        error_message = f"PRG argument provided and file not found: {args.prg}"
        log.error(error_message)
        raise FileNotFoundError(error_message)
    else:
        log.debug("PRG file provided, skipping construction")
        log.debug("Copying PRG file into gram directory")
        shutil.copyfile(args.prg, build_paths['prg'])


## Executes `gram build` backend.
@with_report
def _execute_gramtools_cpp_build(report, action, build_paths, args):

    command = [
        common.gramtools_exec_fpath,
        'build',
        '--gram', build_paths['gram_dir'],
        '--kmer-size', str(args.kmer_size),
        '--max-read-size', str(args.max_read_length),
        '--max-threads', str(args.max_threads),
    ]

    if args.all_kmers:
        command.append('--all-kmers')

    if args.debug:
        command += ['--debug']
    command_str = ' '.join(command)

    log.debug('Executing command:\n\n%s\n', command_str)

    current_working_directory = os.getcwd()
    log.debug('Using current working directory:\n%s',
              current_working_directory)

    process_handle = subprocess.Popen(command_str,
                                      cwd=current_working_directory,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE,
                                      shell=True,
                                      env={'LD_LIBRARY_PATH': common.lib_paths})

    process_result = common.handle_process_result(process_handle)
    command_result, entire_stdout = process_result

    # Add extra reporting
    report[action] = collections.OrderedDict([
        ('command', command_str),
        ('stdout', entire_stdout),
    ])
    if command_result == False:
        raise Exception("Error running gramtools build.")


def _save_report(start_time,
                 report,
                 command_paths,
                 command_hash_paths,
                 report_file_path):
    end_time = str(time.time()).split('.')[0]
    _, report_dict = version.report()
    current_working_directory = os.getcwd()

    _report = collections.OrderedDict([
        ('start_time', start_time),
        ('end_time', end_time),
        ('total_runtime', int(end_time) - int(start_time)),
    ])
    _report.update(report)
    _report.update(collections.OrderedDict([
        ('current_working_directory', current_working_directory),
        ('paths', command_paths),
        ('path_hashes', command_hash_paths),
        ('version_report', report_dict),
    ]))

    with open(report_file_path, 'w') as fhandle:
        json.dump(_report, fhandle, indent=4)


def get_next_valid_record(vcf_reader):
    try:
        next_record = next(vcf_reader, None)
    except Exception:
        next_record = get_next_valid_record(vcf_reader)
    return next_record
