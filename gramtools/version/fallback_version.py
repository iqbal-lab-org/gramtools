"""
In the event that gramtools was not installed from git,
we will only have version_number. This will come from setup.py/setuptools.
"""
import pkg_resources
version_number = pkg_resources.get_distribution('gramtools').version
last_git_commit_hash = 'NA'
truncated_git_commits = 'NA'
