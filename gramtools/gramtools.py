import argparse
import logging

from . import build
from . import kmers
from . import quasimap

try:
    raise ImportError
    from .version import version
except ImportError:
    from .version import fallback_version as version


def setup_logging(level):
    log = logging.getLogger('gramtools')
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(level)
    return log


def _parse_build(common_parser, subparsers):
    build_parser = subparsers.add_parser('build',
                                         parents=[common_parser])
    build_parser.add_argument("--vcf", help="",
                              type=str)
    build_parser.add_argument("--kmer-size", help="",
                              type=int)
    build_parser.add_argument("--reference", help="",
                              type=str)
    build_parser.add_argument("--kmer-region-distance",
                              dest="kmer_region_distance",
                              help="",
                              type=int)
    build_parser.add_argument("--nonvariant-kmers", help="",
                              default=False,
                              action="store_true")
    build_parser.add_argument("--output-fpath", help="",
                              default='',
                              type=str)


def _parse_kmers(common_parser, subparsers):
    kmers_parser = subparsers.add_parser('kmers',
                                         parents=[common_parser])
    kmers_parser.add_argument("--kmer-size", help="",
                              type=int)
    kmers_parser.add_argument("--reference", help="",
                              type=str)
    kmers_parser.add_argument("--kmer-region-distance",
                              dest="kmer_region_distance",
                              help="",
                              type=int)
    kmers_parser.add_argument("--nonvariant-kmers", help="",
                              action="store_true")
    kmers_parser.add_argument("--output-fpath", help="",
                              type=str)


def _parse_quasimap(common_parser, subparsers):
    quasimap_parser = subparsers.add_parser('quasimap',
                                            parents=[common_parser])
    quasimap_parser.add_argument("--gram-files", help="",
                                 type=str)
    quasimap_parser.add_argument("--fastaq", help="",
                                 type=str)
    quasimap_parser.add_argument("--kmer-size", help="",
                                 type=int)


def parse_args():
    root_parser = argparse.ArgumentParser(prog='gramtools')
    root_parser.add_argument("--version", help="",
                             action="store_true")
    subparsers = root_parser.add_subparsers(title='subcommands',
                                            dest='subparser_name',
                                            metavar='{build, kmers, quasimap}')

    common_parser = subparsers.add_parser('common', add_help=False)
    common_parser.add_argument("--debug", help="",
                               action="store_true")
    common_parser.add_argument("--profile", help="",
                               action="store_true")

    _parse_build(common_parser, subparsers)
    _parse_kmers(common_parser, subparsers)
    _parse_quasimap(common_parser, subparsers)

    arguments = root_parser.parse_args()
    return root_parser, arguments


def report_version(log):
    log.info("Latest commit hash:\n%s", version.latest_commit)
    log.info("Current branch: %s", version.current_branch)

    if version.commit_log == 'NA':
        commits = version.commit_log
    else:
        commits = version.commit_log.split('*****')[1:]
        commits = '\n'.join(commits)
    log.info("Truncated commit log:\n%s", commits)


def run():
    root_parser, args = parse_args()

    if hasattr(args, 'debug') and args.debug:
        level = logging.DEBUG
    else:
        level = logging.INFO
    setup_logging(level)
    log = logging.getLogger('gramtools')

    if args.version:
        report_version(log)
        return

    command_switch = {
        'build': build,
        'kmers': kmers,
        'quasimap': quasimap,
    }

    try:
        command = command_switch[args.subparser_name]
    except KeyError:
        root_parser.print_help()
    else:
        command.run(args)


if __name__ == '__main__':
    run()
