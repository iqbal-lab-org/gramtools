all:
	@echo "\
	latest_commit = \"$(shell git rev-parse HEAD)\"\n\
	current_branch = \"$(shell git rev-parse --abbrev-ref HEAD)\"\n\
	commit_log = \"$(shell git log --pretty=format:"*****%h - %an, %ar : %s" -5)\"" > git_version.py;
