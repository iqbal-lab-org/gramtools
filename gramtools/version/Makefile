# We parse setup.py to get the version number and use git for commit status
all:
	# @echo "version_number = $(shell grep 'version=' ../../setup.py | sed 's/version=//' | tr -d ' ,')" > version.py;
	@echo "last_git_commit_hash = \"$(shell git rev-parse HEAD)\"" > commit_version.py;
	@echo "truncated_git_commits = \"$(shell git log --pretty=format:"*****%h - %an, %at : %s" -5)\"" >> commit_version.py;
