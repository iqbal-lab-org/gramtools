"""
Build/load a population reference genome and set it up for quasimapping.
Either a vcf/reference is passed and a prg generated from it, or an existing prg is passed.
Once the prg is stored the back-end `build` routine is called, producing the encoded prg, its fm-index, and other supporting data structures.
"""
import shutil
import logging
import collections

import cluster_vcf_records

from gramtools.commands import common, report
from . import vcf_to_prg_string
from . import command_setup

log = logging.getLogger("gramtools")


def run(args):
    build_paths = command_setup.setup_files(args)

    log.info("Start process: build")
    build_report = report.new_report()

    _prepare_prg(build_report, build_paths, args)
    _execute_gramtools_cpp_build(build_report, "gramtools_build", build_paths, args)

    log.debug("Computing sha256 hash of project paths")
    command_hash_paths = common.hash_command_paths(build_paths)

    build_report.update(collections.OrderedDict({"kmer_size": args.kmer_size}))

    report._save_report(build_report, build_paths, command_hash_paths)
    log.info(f"Success! Build process report in {build_paths.report}")


def _count_vcf_record_lines(vcf_file_path):
    num_recs = 0
    with open(vcf_file_path) as f_in:
        for line in f_in:
            if line[0] != "#":
                num_recs += 1
    return num_recs


def _prepare_prg(build_report, build_paths, args):
    if args.prg is not None:
        _skip_prg_construction(
            build_report, "copy_existing_PRG_string", build_paths, args
        )
    else:
        if args.no_vcf_clustering:
            _skip_cluster_vcf_records(
                build_report, "skip_vcf_record_clustering", build_paths
            )
        else:
            # Update the vcf path to a combined vcf from all those provided.
            # We also do this if only a single one is provided, to deal with overlapping records.
            _cluster_vcf_records(build_report, "vcf_record_clustering", build_paths)
        _execute_command_generate_prg(
            build_report, "vcf_to_PRG_string_conversion", build_paths
        )


@report.with_report
def _skip_cluster_vcf_records(report, action, build_paths):
    if len(build_paths.input_vcfs) > 1:
        log.error(
            "If you ask for no clustering, please provide a single vcf file as input"
        )
        exit(1)
    shutil.copy(build_paths.input_vcfs[0], build_paths.built_vcf)


## Combines records in one or more vcf files using external python utility.
# Records where the REFs overlap are merged together and all possible haplotypes enumerated.
# New path to 'vcf' is path to the combined vcf
@report.with_report
def _cluster_vcf_records(report, action, build_paths):
    log.info(f"Running {action} on {build_paths.input_vcfs}.")

    cluster = cluster_vcf_records.vcf_clusterer.VcfClusterer(
        [str(i) for i in build_paths.input_vcfs],
        str(build_paths.ref),
        str(build_paths.built_vcf),
        max_alleles_per_cluster=5000,
    )
    cluster.run()


## Calls utility that converts a vcf and fasta reference into a linear prg.
@report.with_report
def _execute_command_generate_prg(report, action, build_paths):
    built_vcf = build_paths.built_vcf
    log.info(f"Running {action} on {built_vcf}")

    converter = vcf_to_prg_string.Vcf_to_prg(
        built_vcf, build_paths.ref, build_paths.prg, mode="normal"
    )
    converter._write_bytes()
    converter._write_coordinates()

    num_recs_in_vcf = _count_vcf_record_lines(built_vcf)
    assert num_recs_in_vcf == converter.num_sites, log.error(
        f"Mismatch between number of vcf records in {built_vcf}"
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
    shutil.copyfile(args.prg, build_paths.prg)

    # Write coordinates file
    with open(build_paths.coords_file, "w") as genome_file:
        if args.reference != "None":  # empty file signals no segments
            for rec_id, rec_size in common.load_fasta(
                args.reference, sizes_only=True
            ).items():
                line = f"{rec_id}\t{rec_size}\n"
                genome_file.write(line)


## Executes `gram build` backend.
@report.with_report
def _execute_gramtools_cpp_build(build_report, action, build_paths, args):
    log.info("Running backend build")
    command = [
        common.gramtools_exec_fpath,
        "build",
        "--gram",
        str(build_paths.gram_dir),
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

    command_result, entire_stdout = common.run_subprocess(command)

    # Add extra reporting
    build_report["processes"][action] = collections.OrderedDict(
        [("command", "".join(command)), ("stdout", entire_stdout)]
    )
    if command_result == False:
        raise Exception("Error running gramtools build.")
