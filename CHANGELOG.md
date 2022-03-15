# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
* Support for building a genome graph from multiple MSAs and/or genome sub-graphs built by [make_prg][make_prg].
  Addresses [#130][130].
* CONTRIBUTING.md document

### Changed
* Dependencies: added [make_prg][make_prg] and pybedtools, updated biopython version.

## [1.9.0] - 25/01/2022

### Added

* This Changelog
* Static executable statically-linked against glibc/libgcc/libstd++ [[#101][101],[#163][163]]
* Script to build the gramtools backend (build.sh). Simplifies CI (docker image building), removes
  hacking into setup.py, and simplifies compiling on user machine.

### Changed

* Switched from setup.py to pyproject.toml and setup.cfg [[#164][164]],
  following PEP517

### Removed

* Retired variant-aware kmer indexing code [[#159][159]]
* No longer install `py-cortex-api` by default. If `gramtools discover` is used,
  its installation is requested

### Fixed

* Variant rebasing in `discover` was incompletely tested- simplified algorithm and added tests
* Infinite loop when no reads mapped to genome graph

## [1.8.0] - 18/05/2021

### Added

* Use [conan][conan] dependency manager for available dependencies (boost, gtest, json) [[#158][158]]
* Region-based (chr:start-stop) subgraph visualisation (libgramtools/submods/visualise_prg)

### Changed

* Bump CMake requirements to 3.1.2
* Bump boost libraries to 1.69.0


[unreleased]: https://github.com/iqbal-lab-org/gramtools/compare/1.9.0...HEAD
[1.9.0]: https://github.com/iqbal-lab-org/gramtools/releases/tag/v1.9.0
[1.8.0]: https://github.com/iqbal-lab-org/gramtools/releases/tag/v1.8.0
[101]: https://github.com/iqbal-lab-org/gramtools/issues/101
[130]: https://github.com/iqbal-lab-org/gramtools/issues/130
[158]: https://github.com/iqbal-lab-org/gramtools/issues/158
[159]: https://github.com/iqbal-lab-org/gramtools/issues/159
[163]: https://github.com/iqbal-lab-org/gramtools/issues/163
[164]: https://github.com/iqbal-lab-org/gramtools/issues/164

[conan]: https://conan.io/
[make_prg]: https://github.com/iqbal-lab-org/make_prg
