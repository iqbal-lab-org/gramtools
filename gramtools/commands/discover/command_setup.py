from gramtools import py_cortex_api_message


def setup_parser(common_parser, subparsers):
    parser = subparsers.add_parser(
        "discover",
        parents=[common_parser],
        usage=f"%(prog)s [-h] [options] \n{py_cortex_api_message}",
    )

    parser.add_argument(
        "-i",
        "--genotype_dir",
        help="Input directory."
        "Is an output directory of gramtools `genotype` command.",
        dest="geno_dir",
        type=str,
        required=True,
    )

    parser.add_argument(
        "-o",
        "--discovery_dir",
        help="Directory for outputs,"
        "which are:\n i)the discovered variants expressed against"
        " the inferred personalised reference \n ii)the final variants"
        " expressed against the original genome reference."
        " The latter is useful for augmenting the PRG (= using as input in gramtools `build`).",
        dest="disco_dir",
        type=str,
        required=True,
    )
