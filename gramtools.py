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


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--version", help="",
                        action="store_true")
    parser.add_argument("--build", help="",
                        action="store_true")
    parser.add_argument("--quasimap", help="",
                        action="store_true")

    parser.add_argument("--debug", help="",
                        action="store_true")
    parser.add_argument("--profile", help="",
                        action="store_true")

    parser.add_argument("--gram-files", help="",
                        type=str)
    parser.add_argument("--reference", help="",
                        type=str)
    parser.add_argument("--fastaq", help="",
                        type=str)
    parser.add_argument("--vcf", help="",
                        type=str)

    parser.add_argument("--ksize", help="",
                        type=int)
    parser.add_argument("--kmer-region-distance",
                        dest="kmer_region_distance",
                        help="",
                        type=int)

    args = parser.parse_args()
    return args


def report_version(log):
    log.info("Latest commit hash:\n%s", git_version.latest_commit)
    log.info("Current branch: %s", git_version.current_branch)

    commits = git_version.commit_log.split('*****')[1:]
    log.info("Truncated commit log:\n%s", '\n'.join(commits))


if __name__ == '__main__':
    args = parse_args()

    if args.debug:
        level = logging.DEBUG
    else:
        level = logging.INFO
    setup_logging(level)
    log = logging.getLogger('gramtools')

    if args.version:
        report_version(log)

    elif args.build:
        build.run(args)

    elif args.quasimap:
        quasimap.run(args)
