import json
import logging
import argparse
import collections

from . import build
from . import kmers
from . import simulate
from . import quasimap

try:
    from .version import version
except ImportError:
    from .version import fallback_version as version


def _setup_logging(level):
    log = logging.getLogger('gramtools')
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(level)
    return log


root_parser = argparse.ArgumentParser(prog='gramtools')
commands = collections.OrderedDict([
    ('build', build),
    ('kmers', kmers),
    ('simulate', simulate),
    ('quasimap', quasimap),
])


def _parse_args():
    root_parser.add_argument('--version', help='',
                             action='store_true')
    metavar = '{{{commands}}}'.format(commands=', '.join(commands.keys()))
    subparsers = root_parser.add_subparsers(title='subcommands',
                                            dest='subparser_name',
                                            metavar=metavar)

    common_parser = subparsers.add_parser('common', add_help=False)
    common_parser.add_argument('--debug', help='',
                               action='store_true')
    common_parser.add_argument('--profile', help='',
                               action='store_true')

    for command in commands.values():
        command.parse_args(common_parser, subparsers)

    arguments = root_parser.parse_args()
    return arguments


def _report_version(log):
    if version.truncated_git_commits == 'NA':
        commits = []
    else:
        commits = version.truncated_git_commits.split('*****')[1:]
        commits = [x.strip() for x in commits]

    report = collections.OrderedDict([
        ('version_number', version.version_number),
        ('last_git_commit_hash', version.last_git_commit_hash),
        ('current_git_branch', version.current_git_branch),
        ('truncated_git_commits', commits),
    ])
    report_json = json.dumps(report, indent=4)
    print(report_json)


def _get_log(args):
    if hasattr(args, 'debug') and args.debug:
        level = logging.DEBUG
    else:
        level = logging.INFO
    _setup_logging(level)
    log = logging.getLogger('gramtools')
    return log


def run():
    args = _parse_args()
    log = _get_log(args)

    if args.version:
        _report_version(log)
        return

    try:
        command = commands[args.subparser_name]
    except KeyError:
        root_parser.print_help()
    else:
        command.run(args)


if __name__ == '__main__':
    run()
