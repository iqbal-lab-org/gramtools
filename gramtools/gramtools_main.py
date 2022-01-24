## @file
# The entry point for the `gramtools` software.
import logging
import argparse
import collections
from subprocess import run as sp_run
from pathlib import Path
import sys

from gramtools import version, gramtools_exec_fpath
from gramtools.commands.build import command_setup as build_setup, build
from gramtools.commands.genotype import command_setup as genotype_setup, genotype
from gramtools.commands.discover import command_setup as discovery_setup, discover
from gramtools.commands.simulate import simulate


def _setup_logging(args):
    log = logging.getLogger("gramtools")
    log.propagate = False  # Do not pass to ancestor loggers
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(asctime)s %(name)-12s %(levelname)-8s %(message)s")
    handler.setFormatter(formatter)
    log.addHandler(handler)

    if hasattr(args, "debug") and args.debug:
        level = logging.DEBUG
    else:
        level = logging.INFO
    log.setLevel(level)


root_parser = argparse.ArgumentParser(prog="gramtools")
command_setups = collections.OrderedDict(
    [
        ("build", build_setup),
        ("genotype", genotype_setup),
        ("discover", discovery_setup),
    ]
)
commands = collections.OrderedDict(
    [
        ("build", build),
        ("genotype", genotype),
        ("discover", discover),
        ("simulate", simulate),
    ]
)


def _setup_parser():
    root_parser.add_argument("--version", help="", action="store_true")
    metavar = "{{{commands}}}".format(commands=", ".join(commands.keys()))
    subparsers = root_parser.add_subparsers(
        title="subcommands", dest="subparser_name", metavar=metavar
    )

    # Will add a --debug mode for all commands
    common_parser = subparsers.add_parser("common", add_help=False)
    common_parser.add_argument(
        "--debug", help="Verbose logging of actions taken", action="store_true"
    )
    common_parser.add_argument(
        "--force",
        help="Erase already existing output directory/files",
        action="store_true",
    )

    for command_setup in command_setups.values():
        command_setup.setup_parser(common_parser, subparsers)
    simulate.setup_parser(common_parser, subparsers)


def check_gram_binary():
    if not Path(gramtools_exec_fpath).exists():
        print(
            f"gramtools backend expected at: {gramtools_exec_fpath} but not found.",
            file=sys.stderr,
        )
        exit(1)
    process = sp_run(
        [gramtools_exec_fpath], capture_output=True, universal_newlines=True
    )
    if process.returncode != 0:
        print(
            f"The gramtools backend at {gramtools_exec_fpath} does not seem to work.\n\n"
            f"Stdout: \n{process.stdout.strip()}\n\n",
            f"Stderr: \n{process.stderr.strip()}\n",
            file=sys.stderr,
        )
        exit(1)


def run():
    _setup_parser()
    args = root_parser.parse_args()

    check_gram_binary()
    _setup_logging(args)
    if args.version:
        report_json, _ = version.report()
        print(report_json)
        return

    if args.subparser_name is None:
        root_parser.print_help()
        exit(1)
    command = commands[args.subparser_name]
    command.run(args)


if __name__ == "__main__":
    run()
