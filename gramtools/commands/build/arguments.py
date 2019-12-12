import argparse


def setup_build_parser(common_parser, subparsers):
    parser = subparsers.add_parser("build", parents=[common_parser])
    parser.add_argument(
        "--gram-dir",
        "--gram-directory",
        help="",
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
        help="A prg can be passed in directly instead of a vcf/reference combination.",
        type=str,
    )

    parser.add_argument(
        "--kmer-size",
        help="Kmer size for indexing the prg. Defaults to 10. "
        "Higher k speeds quasimapping. Currently capped at 14 (268 million kmers), because of all kmers enumeration.",
        type=int,
        default=10,
        required=False,
    )

    parser.add_argument(
        "--max-read-length",
        help=argparse.SUPPRESS,
        # help="Used to determine which kmers overlap variant sites",
        type=int,
        default=150,
        required=False,
    )

    parser.add_argument("--max-threads", help="", type=int, default=1, required=False)


def _check_build_args(args):
    no_prg = args.prg is None
    no_vcf_and_no_ref = args.reference is None and args.vcf is None
    not_both_vcf_and_ref = args.reference is None or args.vcf is None

    if no_prg and no_vcf_and_no_ref:
        print(
            "Please provide known variation through either: \n* --prg \n* --vcf and --reference"
        )
        exit(1)

    if args.kmer_size > 14:
        print(
            "--kmer-size must be 14 or less, because indexing currently produces all kmers of given size."
        )
        exit(1)

    if not no_prg:
        return

    if not_both_vcf_and_ref:
        print("Please provide both --reference and --vcf")
        exit(1)
