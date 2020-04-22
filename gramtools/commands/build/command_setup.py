import argparse
from ..paths import ProjectPaths, BuildPaths
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
        "--ref",
        "--reference",
        help="Reference genome. Used to build non-variant parts of the prg (--vcf) "
        "or to book-keep chromosome IDs and coordinates (--prg).",
        type=str,
        dest="reference",
        required=True,
    )

    variation = parser.add_mutually_exclusive_group(required=True)
    variation.add_argument(
        "--vcf",
        help="File(s) containing variant information to capture in the prg.",
        nargs="+",
        action="append",
        type=str,
    )

    variation.add_argument(
        "--prg",
        help="A prg can be passed in directly instead of a vcf. "
        "Use this for prgs with nesting.",
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
    parser.add_argument(
        "--max_threads", help=argparse.SUPPRESS, type=int, default=1, required=False
    )

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


def setup_files(args) -> BuildPaths:
    """
    We also do some extra argument checking here.
    """
    build_paths = BuildPaths(args.gram_dir, args.force)
    build_paths.setup()

    if args.kmer_size > 14:
        err_message = "--kmer-size must be 14 or less, because indexing currently produces all kmers of given size."
        build_paths.raise_error(err_message)

    if args.vcf is not None:
        build_paths.ready_ref_and_vcf(args.reference, args.vcf)

    return build_paths
