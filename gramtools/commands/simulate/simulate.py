import logging
import time

from gramtools import gramtools_exec_fpath
from gramtools.commands.paths import SimulatePaths
from gramtools.commands import common

log = logging.getLogger("gramtools")


def setup_parser(common_parser, subparsers):
    parser = subparsers.add_parser("simulate", parents=[common_parser])
    parser.add_argument(
        "--prg",
        help="A prg as made by (or passed to) gramtools build",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-n",
        "--max_num_paths",
        help="Number of paths through the prg to simulate. \n"
        "Duplicates are removed, so this is an upper bound.",
        type=int,
        required=False,
        default=100,
    )
    parser.add_argument(
        "--sample_id",
        help="A name for your sampled paths.\n"
        "Prefixes the output filenames and sample IDs in output files.",
        required=False,
        default="sim",
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        help="directory containing outputs",
        type=str,
        required=False,
        default=".",
    )
    parser.add_argument(
        "-i",
        "--induce_genotypes",
        help="Input multifasta to produce genotypes of.\n"
        "Fails if any input sequence is not found in the graph.\n "
        "If multiple paths correspond to one sequence (ie graph is ambiguous), "
        "uses one path.\n"
        "Paths can not fully consume the sequence; if there are several, the longest-consuming"
        " path is used.",
        type=str,
        required=False,
        default="",
    )


def run(args):
    simu_paths = SimulatePaths(
        args.output_dir, args.sample_id, args.prg, args.induce_genotypes, args.force
    )
    simu_paths.setup()

    log.info("Start process: simulate")
    start_time = str(time.time()).split(".")[0]

    _execute_command_cpp_simulate(simu_paths, args)

    log.info("End process: simulate")


def _execute_command_cpp_simulate(simu_paths, args):
    input_multifasta = list()
    if hasattr(simu_paths, "input_multifasta"):
        input_multifasta.extend(["--i", str(simu_paths.input_multifasta)])
    command = [
        gramtools_exec_fpath,
        "simulate",
        "--prg",
        str(simu_paths.prg_fpath),
        "--n",
        str(args.max_num_paths),
        "--sample_id",
        args.sample_id,
        "--o",
        str(simu_paths.output_dir),
    ] + input_multifasta

    if args.debug:
        command += ["--debug"]

    command_result = common.run_subprocess(command)

    if not command_result.success:
        raise Exception(f"Error running gramtools simulate:\n{command_result.stderr}")
