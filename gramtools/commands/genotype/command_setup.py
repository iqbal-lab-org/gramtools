from ..paths import ProjectPaths, GenotypePaths


def setup_parser(common_parser, subparsers):
    parser = subparsers.add_parser("genotype", parents=[common_parser])
    parser.add_argument(
        "-i",
        "--gram_dir",
        help="Directory containing outputs from gramtools `build`",
        dest="gram_dir",
        type=str,
        required=True,
    )

    parser.add_argument(
        "-o",
        "--genotype_dir",
        help="Directory to hold this command's outputs.",
        type=str,
        dest="geno_dir",
        required=True,
    )

    parser.add_argument(
        "--reads",
        help="One or more read files.\n"
        "Valid formats: fastq, sam/bam/cram, fasta, txt; compressed or uncompressed; fuzzy extensions (eg fq, fsq for fastq).\n"
        "Read files can be given after one or several '--reads' argument:"
        "Eg '--reads rf_1.fq rf_2.fq.gz --reads rf_3.bam '",
        nargs="+",
        action="append",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--ploidy",
        help="The expected ploidy of the sample.\n" "Default: diploid",
        choices=["haploid", "diploid"],
        required=False,
        default="diploid",
    )

    parser.add_argument(
        "--max_threads",
        help="Run with more threads than the default of one.",
        type=int,
        default=1,
        required=False,
    )

    parser.add_argument(
        "--seed",
        help="Fixing the seed will produce the same read mappings across different runs."
        "By default, seed is randomly generated so this is not the case.",
        type=int,
        default=0,
        required=False,
    )

    parser.add_argument(
        "--force", help="Erase pre-existing output directory", action="store_true"
    )


## Make quasimap-related file and directory paths.
def setup_files(args) -> ProjectPaths:
    geno_paths = GenotypePaths(args.geno_dir, args.force)
    geno_paths.initial_setup()

    geno_paths._link_to_build(args.gram_dir)
    geno_paths._link_to_reads(args.reads)

    return geno_paths
