all:
	@echo "\
	const std::string latest_commit = \"$(shell git rev-parse HEAD)\";\n\
	const std::string current_branch = \"$(shell git rev-parse --abbrev-ref HEAD)\";\n\
	const std::string commit_log = \"$(shell git log --pretty=format:"*****%h - %an, %ar : %s" -5)\";" > git_version.hpp;
