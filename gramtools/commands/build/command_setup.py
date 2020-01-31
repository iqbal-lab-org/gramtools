import argparse
from ..paths import ProjectPaths, BuildPaths
from pathlib import Path
import logging

log = logging.getLogger("gramtools")


def setup_parser(common_parser, subparsers):
    parser = subparsers.add_parser("build", parents=[common_parser])
    parser.add_argument(
        "-o",
        "--gram_dir",
        help="Directory containing the built prg.",
        dest="gram_dir",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--vcf",
        help="File containing variant information to capture in the prg.",
        nargs="+",
        action="append",
        type=str,
    )

    parser.add_argument(
        "--reference",
        help="Reference the vcf file refers to, used to build non-variant parts of the prg.",
        type=str,
    )
    parser.add_argument(
        "--prg",
        help="A prg can be passed in directly instead of a combination of vcf/reference.",
        type=str,
    )

    parser.add_argument(
        "--kmer_size",
        help="Kmer size for indexing the prg. Defaults to 10. "
        "Higher k speeds quasimapping. Currently capped at 14 (268 million kmers), because of all kmers enumeration.",
        type=int,
        default=10,
        required=False,
    )

    # TODO: there is no multi-threading yet
    # parser.add_argument("--max_threads", help="", type=int, default=1, required=False)

    # Hidden arguments, for legacy/special uses (minos)
    parser.add_argument(
        "--max_read_length",
        help=argparse.SUPPRESS,
        # help="Used to determine which kmers overlap variant sites",
        type=int,
        default=150,
        required=False,
    )

    parser.add_argument(
        "--no_vcf_clustering",
        help=argparse.SUPPRESS,
        # help="Do not run vcf clustering on the input vcfs",
        action="store_true",
    )


def setup_files(args) -> ProjectPaths:
    """
    We also do some extra argument checking here.
    """
    build_paths = BuildPaths(args.gram_dir)

    no_prg = args.prg is None
    no_vcf_and_no_ref = args.reference is None and args.vcf is None
    not_both_vcf_and_ref = args.reference is None or args.vcf is None

    if no_prg and no_vcf_and_no_ref:
        err_message = "Please provide known variation through either: \n* --prg \n* --vcf and --reference"
        build_paths.raise_error(err_message)

    if args.kmer_size > 14:
        err_message = "--kmer-size must be 14 or less, because indexing currently produces all kmers of given size."
        build_paths.raise_error(err_message)

    if not no_prg:
        build_paths.check_exists(Path(args.prg))
        if not no_vcf_and_no_ref:
            log.warning(
                "You have provided a --prg and --reference/--vcf. Building using the --prg argument only."
            )
        return build_paths

    if not_both_vcf_and_ref:
        err_message = "Please provide both --reference and --vcf"
        build_paths.raise_error(err_message)

    # We have been given ref and vcf
    build_paths.ready_ref_and_vcf(args.reference, args.vcf)

    return build_paths
