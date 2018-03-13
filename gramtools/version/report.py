import json
import collections

try:
    from . import version
    if not hasattr(version, 'last_git_commit_hash'):
        raise ImportError
except ImportError:
    from . import fallback_version as version


def report():
    if version.truncated_git_commits == 'NA':
        commits = []
    else:
        commits = version.truncated_git_commits.split('*****')[1:]
        commits = [x.strip() for x in commits]

    report_dict = collections.OrderedDict([
        ('version_number', version.version_number),
        ('last_git_commit_hash', version.last_git_commit_hash),
        ('truncated_git_commits', commits),
    ])
    report_json = json.dumps(report_dict, indent=4)
    return report_json, report_dict
