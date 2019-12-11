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
        help="Kmer size for indexing the prg. Defaults to 5.",
        type=int,
        default=5,
        required=False,
    )

    # The current default behaviour is to extract only relevant kmers from prg.
    parser.add_argument(
        "--all-kmers",
        help="Whether or not all kmers of given size should be indexed.\n"
        "When this flag is not used, only kmers overlapping variant sites in prg will be indexed.",
        action="store_true",
        required=False,
    )

    parser.add_argument(
        "--max-read-length",
        help="Used to determine which kmers overlap variant sites. Only needed if --all-kmers flag is off.",
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
        log.error(
            "Please provide genetic variation via either --prg or --vcf/--reference"
        )
        exit(1)
    if not no_prg:
        return
    if not_both_vcf_and_ref:
        log.error("Please provide both --reference and --vcf")
        exit(1)
