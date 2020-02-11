## @file
# Executes `gram genotype` backend which:
# - maps reads to prg using vBWT-based (quasi)mapping and annotates it with coverage
# - genotypes the prg using the annotated coverage
import os
import time
import json
import logging
import subprocess
import collections

from ... import common
from .. import report
from . import command_setup
from ..paths import GenotypePaths

log = logging.getLogger("gramtools")


def run(args):
    log.info("Start process: genotype")

    start_time = str(time.time()).split(".")[0]
    geno_paths = command_setup.setup_files(args)

    build_report = _load_build_report(geno_paths)

    kmer_size = build_report["kmer_size"]
    setattr(args, "kmer_size", kmer_size)

    geno_report = collections.OrderedDict({"processes": collections.OrderedDict()})
    _execute_command(geno_report["processes"], "gramtools_genotype", geno_paths, args)

    _check_read_stats(geno_report["processes"], "check_read_stats", geno_paths)

    log.debug("Computing sha256 hash of project paths")
    command_hash_paths = common.hash_command_paths(geno_paths)

    log.debug("Saving command report:\n%s", geno_paths.report)
    success_status = {"success": geno_report["processes"].pop("success")}
    geno_report = {**success_status, **geno_report}
    report._save_report(
        start_time, geno_report, geno_paths, command_hash_paths, geno_paths.report
    )

    if geno_report["success"] is False:
        log.error(f"Unsuccessful genotyping. Process reported to {geno_paths.report}")
        exit(1)
    else:
        log.info(f"Success! Genotyping process report in {geno_paths.report}")


def _load_build_report(geno_paths):
    build_path = geno_paths.gram_dir / "build_report.json"

    if not build_path.exists():
        log.error(
            f"Build report not found: {build_path}. Try re-running gramtools `build`?"
        )
        exit(1)

    with open(build_path) as fhandle:
        build_report = json.load(fhandle)
        if not build_report["success"]:
            log.error(f"Build was not completed successfully: see {build_path}")
            exit(1)

    return build_report


@report.with_report
def _execute_command(geno_report, action, geno_paths, args):

    command = [
        common.gramtools_exec_fpath,
        "genotype",
        "--gram_dir",
        str(geno_paths.gram_dir),
        "--reads",
        "".join([str(read_file) for read_file in geno_paths.reads_files]),
        "--ploidy",
        args.ploidy,
        "--kmer_size",
        str(args.kmer_size),
        "--genotype_dir",
        str(geno_paths.geno_dir),
        "--max_threads",
        str(args.max_threads),
        "--seed",
        str(args.seed),
    ]

    command_str = " ".join(command)
    log.debug("Executing command:\n\n%s\n", command_str)

    current_working_directory = os.getcwd()
    process_handle = subprocess.Popen(
        command_str,
        cwd=current_working_directory,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True,
        env={"LD_LIBRARY_PATH": common.lib_paths},
    )

    command_result, entire_stdout = common.handle_process_result(process_handle)
    log.debug("Output run directory:\n%s", geno_paths.geno_dir)

    # Add extra reporting
    geno_report[action] = collections.OrderedDict(
        [("command", command_str), ("stdout", entire_stdout)]
    )
    if command_result == False:
        raise Exception("Error running gramtools genotype.")


@report.with_report
def _check_read_stats(geno_report, action, geno_paths: GenotypePaths):
    """
    Get the read statistics; report if most variant sites have no coverage.
    """
    with open(geno_paths.read_stats) as f:
        read_stats = json.load(f)

    num_sites_noCov, num_sites_total = (
        read_stats["Read_depth"]["num_sites_noCov"],
        read_stats["Read_depth"]["num_sites_total"],
    )
    if num_sites_noCov / num_sites_total > 0.5:
        log.warning(
            "More than 50% of all variant sites have no coverage ({} out of {})."
            "Possible reasons include: reads not quality-trimmed; low sequencing depth.".format(
                num_sites_noCov, num_sites_total
            )
        )
