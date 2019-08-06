import json
import collections
from . import package_version


def report():
    try:
        from . import commit_version
        commits = commit_version.truncated_git_commits.split('*****')[1:]
        commits = [x.strip() for x in commits]
        last_commit_hash = commit_version.last_git_commit_hash
    except ImportError:
        commits = []
        last_commit_hash = "NA"

    report_dict = collections.OrderedDict([
        ('version_number', package_version.__version__),
        ('last_git_commit_hash', last_commit_hash),
        ('truncated_git_commits', commits),
    ])
    report_json = json.dumps(report_dict, indent=4)
    return report_json, report_dict
