import argparse
import logging

from py_interface import build, quasimap
from py_interface.git_version import git_version


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
                                            metavar='{build,quasimap}')

    common_parser = subparsers.add_parser('common', add_help=False)
    common_parser.add_argument("--debug", help="",
                               action="store_true")
    common_parser.add_argument("--profile", help="",
                               action="store_true")

    _parse_build(common_parser, subparsers)
    _parse_quasimap(common_parser, subparsers)

    arguments = root_parser.parse_args()
    return arguments


def report_version(log):
    log.info("Latest commit hash:\n%s", git_version.latest_commit)
    log.info("Current branch: %s", git_version.current_branch)

    commits = git_version.commit_log.split('*****')[1:]
    log.info("Truncated commit log:\n%s", '\n'.join(commits))


if __name__ == '__main__':
    args = parse_args()

    if hasattr(args, 'debug') and args.debug:
        level = logging.DEBUG
    else:
        level = logging.INFO
    setup_logging(level)
    log = logging.getLogger('gramtools')

    if args.version:
        report_version(log)

    elif args.subparser_name == 'build':
        build.run(args)

    elif args.subparser_name == 'quasimap':
        quasimap.run(args)
