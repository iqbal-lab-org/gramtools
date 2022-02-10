"""
Build/load a population reference genome (prg) and set it up for quasimapping.
Either a vcf/reference is passed and a prg generated from it, or an existing prg is passed.
Once the prg is stored the back-end `build` routine is called, producing the encoded prg, its fm-index, and other supporting data structures.
"""
import shutil
import logging
import collections

from gramtools import gramtools_exec_fpath
from gramtools.commands import common, report
from gramtools.commands.paths import BuildPaths
from . import command_setup

log = logging.getLogger("gramtools")


def run(args):
    build_paths = command_setup.setup_files(args)

    log.info("Start process: build")
    build_report = report.new_report()

    prepare_prg(build_report, build_paths, args)
    execute_gramtools_cpp_build(build_report, "gramtools_build", build_paths, args)

    log.debug("Computing sha256 hash of project paths")
    command_hash_paths = common.hash_command_paths(build_paths)

    build_report.update(collections.OrderedDict({"kmer_size": args.kmer_size}))

    report._save_report(build_report, build_paths, command_hash_paths)
    log.info(f"Success! Build process report in {build_paths.report}")


def prepare_prg(build_report, build_paths, args):
    chrom_seqs = common.load_fasta(args.reference)
    common.write_coordinates_file(chrom_seqs, build_paths.coords_file)

    if args.prg is not None:
        skip_prg_construction(
            build_report, "copy_existing_PRG_string", build_paths, args
        )
    else:
        build_from_vcfs(build_report, "vcf_to_PRG_string_conversion", build_paths, args)


@report.with_report
def skip_prg_construction(report, action, build_paths, args):
    """Checks prg file exists and copies it to gram directory."""

    log.debug("PRG file provided, skipping construction")
    log.debug("Copying PRG file into gram directory")
    shutil.copyfile(args.prg, build_paths.prg)


@report.with_report
def execute_gramtools_cpp_build(build_report, action, build_paths, args):
    """Executes `gram build` backend."""

    log.info("Running backend build")
    command = [
        gramtools_exec_fpath,
        "build",
        "--gram_dir",
        str(args.gram_dir),
        "--ref",
        str(args.reference),
        "--kmer_size",
        str(args.kmer_size),
        "--max_threads",
        str(args.max_threads),
        "--all_kmers",  # Currently always build all kmers of given size
    ]

    if args.debug:
        command += ["--debug"]

    command_result = common.run_subprocess(command)

    # Add extra reporting
    build_report["processes"][action] = collections.OrderedDict(
        [
            ("command", " ".join(command)),
            ("stdout", command_result.stdout.splitlines()),
            ("stderr", command_result.stderr.splitlines()),
        ]
    )
    if not command_result.success:
        raise Exception(f"while running backend build:\n{command_result.stderr}")
