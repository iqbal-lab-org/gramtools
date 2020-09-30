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
        " eg '--reads rf_1.fq rf_2.fq.gz --reads rf_3.bam '",
        nargs="+",
        action="append",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--sample_id",
        help="A name for your dataset.\n" "Appears in the genotyping outputs.",
        required=True,
    )

    parser.add_argument(
        "--ploidy",
        help="The expected ploidy of the sample.\n" "Default: haploid",
        choices=["haploid", "diploid"],
        required=False,
        default="haploid",
    )

    parser.add_argument(
        "--max_threads",
        help="Max number of threads to use. Default: 1.",
        type=int,
        default=1,
        required=False,
    )

    parser.add_argument(
        "--seed",
        help="Fix the seed to produce the same read mappings across different runs."
        " Default: None (seed gets randomly generated).",
        type=int,
        required=False,
    )
