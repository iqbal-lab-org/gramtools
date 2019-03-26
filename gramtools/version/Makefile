all:
	@echo "version_number = \"2.0.0\"" > version.py;
	@echo "last_git_commit_hash = \"$(shell git rev-parse HEAD)\"" >> version.py;
	@echo "truncated_git_commits = \"$(shell git log --pretty=format:"*****%h - %an, %at : %s" -5)\"" >> version.py;
