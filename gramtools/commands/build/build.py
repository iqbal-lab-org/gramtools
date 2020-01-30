## @file
#  Build/load a population reference genome and set it up for quasimapping.
# Either a vcf/reference is passed and a prg generated from it, or an existing prg is passed.
# Once the prg is stored the back-end `build` routine is called, producing the encoded prg, its fm-index, and other supporting data structures.
import os
import time
import shutil
import logging
import subprocess
import collections

import cluster_vcf_records

from ... import common
from ... import paths
from . import vcf_to_prg_string

from . import arguments
from . import report

log = logging.getLogger("gramtools")


def run(args):
    arguments._check_build_args(args)
    log.info("Start process: build")
    start_time = str(time.time()).split(".")[0]

    command_paths = paths.generate_build_paths(args)
    build_report = collections.OrderedDict()

    if args.prg is not None:
        _skip_prg_construction(
            build_report, "copy_existing_PRG_string", command_paths, args
        )
    else:
        if args.no_vcf_clustering:
            command_paths["vcf"] = _skip_cluster_vcf_records(
                build_report, "skip_vcf_record_clustering", command_paths
            )
        else:
            # Update the vcf path to a combined vcf from all those provided.
            # We also do this if only a single one is provided, to deal with overlapping records.
            command_paths["vcf"] = _cluster_vcf_records(
                build_report, "vcf_record_clustering", command_paths
            )
        _execute_command_generate_prg(
            build_report, "vcf_to_PRG_string_conversion", command_paths
        )

    _execute_gramtools_cpp_build(build_report, "gramtools_build", command_paths, args)

    log.debug("Computing sha256 hash of project paths")
    command_hash_paths = common.hash_command_paths(command_paths)

    log.debug("Saving command report:\n%s", command_paths["build_report"])
    _report = collections.OrderedDict(
        [("kmer_size", args.kmer_size), ("max_read_length", args.max_read_length)]
    )
    build_report.update(_report)

    report_file_path = command_paths["build_report"]
    report._save_report(
        start_time, build_report, command_paths, command_hash_paths, report_file_path
    )

    if build_report["success"] is False:
        log.error(
            f'Unsuccessful build. Process reported to {command_paths["build_report"]}'
        )
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


def setup_command_parser(common_parser, subparsers):
    arguments.setup_build_parser(common_parser, subparsers)


@report.with_report
def _skip_cluster_vcf_records(report, action, command_paths):
    if len(command_paths["vcf"]) > 1:
        log.error(
            "If you ask for no clustering, please provide a single vcf file as input"
        )
        exit(1)
    shutil.copy(command_paths["vcf"][0], command_paths["built_vcf"])
    return command_paths["built_vcf"]


## Combines records in one or more vcf files using external python utility.
# Records where the REFs overlap are merged together and all possible haplotypes enumerated.
# New path to 'vcf' is path to the combined vcf
@report.with_report
def _cluster_vcf_records(report, action, command_paths):
    vcf_files = command_paths["vcf"]
    final_vcf_name = command_paths["built_vcf"]
    log.info(f"Running {action} on {str(vcf_files)}.")

    cluster = cluster_vcf_records.vcf_clusterer.VcfClusterer(
        vcf_files,
        command_paths["original_reference"],
        final_vcf_name,
        max_alleles_per_cluster=5000,
    )
    cluster.run()

    return final_vcf_name


## Calls utility that converts a vcf and fasta reference into a linear prg.
@report.with_report
def _execute_command_generate_prg(report, action, build_paths):
    vcf_in = build_paths["vcf"]
    log.info(f"Running {action} on {vcf_in}")

    converter = vcf_to_prg_string.Vcf_to_prg(
        vcf_in,
        build_paths["original_reference"],
        build_paths["prg_string"],
        mode="normal",
    )
    # print(converter.prg_vector)
    converter._write_bytes()

    ## The converter does not produce a vcf
    # Thus the input vcf needs to be 'clean': we need each vcf record to be converted
    # to a variant site in the prg so that later on the genotyping process, which uses vcf, is possible.
    num_recs_in_vcf = _count_vcf_record_lines(vcf_in)
    assert num_recs_in_vcf == converter.num_sites, log.error(
        f"Mismatch between number of vcf records in {vcf_in}"
        f"({num_recs_in_vcf} and number of variant sites in"
        f"PRG string ({converter.num_sites}.\n"
        f"Possible source of error: vcf record clustering does not"
        f"produce non-overlapping records, or conversion utility"
        f" {vcf_to_prg_string.__file__} is broken."
    )


## Checks prg file exists and copies it to gram directory.
@report.with_report
def _skip_prg_construction(report, action, build_paths, args):
    log.debug("PRG file provided, skipping construction")
    log.debug("Copying PRG file into gram directory")
    shutil.copyfile(args.prg, build_paths["prg"])


## Executes `gram build` backend.
@report.with_report
def _execute_gramtools_cpp_build(report, action, build_paths, args):

    command = [
        common.gramtools_exec_fpath,
        "build",
        "--gram",
        build_paths["gram_dir"],
        "--kmer-size",
        str(args.kmer_size),
        "--max-threads",
        str(args.max_threads),
        "--all-kmers",  # Currently always build all kmers of given size
        "--max-read-size",  # TODO: currently only used when --all-kmers is not passed. Maybe retire entirely.
        "150",
    ]

    if args.debug:
        command += ["--debug"]
    command_str = " ".join(command)

    log.debug("Executing command:\n\n%s\n", command_str)

    current_working_directory = os.getcwd()
    log.debug("Using current working directory:\n%s", current_working_directory)

    process_handle = subprocess.Popen(
        command_str,
        cwd=current_working_directory,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True,
        env={"LD_LIBRARY_PATH": common.lib_paths},
    )

    process_result = common.handle_process_result(process_handle)
    command_result, entire_stdout = process_result

    # Add extra reporting
    report[action] = collections.OrderedDict(
        [("command", command_str), ("stdout", entire_stdout)]
    )
    if command_result == False:
        raise Exception("Error running gramtools build.")
