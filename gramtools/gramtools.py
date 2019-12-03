## @file
# The entry point for the `gramtools` software.
import logging
import argparse
import collections

from . import version
from .commands.build import build
from .commands import quasimap
from .commands.infer import infer
from .commands import discover
from .commands import simulate


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
commands = collections.OrderedDict(
    [("build", build), ("quasimap", quasimap), ("infer", infer), ("discover", discover)]
)


def _setup_parser():
    root_parser.add_argument("--version", help="", action="store_true")
    metavar = "{{{commands}}}".format(commands=", ".join(commands.keys()))
    subparsers = root_parser.add_subparsers(
        title="subcommands", dest="subparser_name", metavar=metavar
    )

    common_parser = subparsers.add_parser("common", add_help=False)
    common_parser.add_argument("--debug", help="", action="store_true")

    for command in commands.values():
        command.setup_command_parser(common_parser, subparsers)


def run():
    _setup_parser()
    args = root_parser.parse_args()

    _setup_logging(args)
    if args.version:
        report_json, _ = version.report()
        print(report_json)
        return

    try:
        command = commands[args.subparser_name]
    except KeyError:
        root_parser.print_help()
    else:
        command.run(args)


if __name__ == "__main__":
    run()
