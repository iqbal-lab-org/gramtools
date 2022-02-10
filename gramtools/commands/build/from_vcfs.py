import logging

from cluster_vcf_records.vcf_clusterer import VcfClusterer

from gramtools.commands import report
from gramtools.commands.paths import BuildPaths
from gramtools.commands.build.vcf_to_prg_string import Vcf_to_prg

log = logging.getLogger("gramtools")


def _count_vcf_record_lines(vcf_file_path):
    num_recs = 0
    with open(vcf_file_path) as f_in:
        for line in f_in:
            if line[0] != "#":
                num_recs += 1
    return num_recs


@report.with_report
def _skip_cluster_vcf_records(report, action, build_paths):
    if len(build_paths.input_vcfs) > 1:
        log.error(
            "If you ask for no clustering, please provide a single vcf file as input"
        )
        exit(1)
    shutil.copy(build_paths.input_vcfs[0], build_paths.built_vcf)


@report.with_report
def cluster_vcf_records(report, action, build_paths: BuildPaths):
    """
    Combines records in one or more vcf files using external python utility.
    Records where the REFs overlap are merged together and all possible haplotypes enumerated.
    New path to 'vcf' is path to the combined vcf
    """
    log.info(f"Running {action} on {build_paths.input_vcfs}.")

    cluster = VcfClusterer(
        [str(i) for i in build_paths.input_vcfs],
        str(build_paths.ref),
        str(build_paths.built_vcf),
        max_alleles_per_cluster=5000,
    )
    cluster.run()


@report.with_report
def build_from_vcfs(report, action, build_paths, args):
    """Calls utility that converts a vcf and fasta reference into a linear prg."""
    if args.no_vcf_clustering:
        _skip_cluster_vcf_records(report, "skip_vcf_record_clustering", build_paths)
    else:
        # We also do this if only a single one is provided, to deal with overlapping records.
        cluster_vcf_records(report, "vcf_record_clustering", build_paths)

    built_vcf = build_paths.built_vcf
    log.info(f"Running {action} on {built_vcf}")

    converter = Vcf_to_prg(built_vcf, build_paths.ref, build_paths.prg, mode="normal")
    converter._write_bytes()

    num_recs_in_vcf = _count_vcf_record_lines(built_vcf)
    assert num_recs_in_vcf == converter.num_sites, log.error(
        f"Mismatch between number of vcf records in {built_vcf}"
        f"({num_recs_in_vcf} and number of variant sites in"
        f"PRG string ({converter.num_sites}.\n"
        f"Please report this to developers."
    )
