all:
	@echo "latest_commit = \"$(shell git rev-parse HEAD)\"" > version.py;
	@echo "current_branch = \"$(shell git rev-parse --abbrev-ref HEAD)\"" >> version.py;
	@echo "commit_log = \"$(shell git log --pretty=format:"*****%h - %an, %ar : %s" -5)\"" >> version.py;
