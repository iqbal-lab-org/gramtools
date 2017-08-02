all:
	@echo "const std::string latest_commit = \"$(shell git rev-parse HEAD)\";" > git_version.hpp;
	@echo "const std::string current_branch = \"$(shell git rev-parse --abbrev-ref HEAD)\";" >> git_version.hpp;
	@echo "const std::string commit_log = \"$(shell git log --pretty=format:"*****%h - %an, %ar : %s" -5)\";" >> git_version.hpp;
