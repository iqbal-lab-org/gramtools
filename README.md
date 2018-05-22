[![Build Status](https://travis-ci.org/iqbal-lab-org/gramtools.svg?branch=master)](https://travis-ci.org/iqbal-lab-org/gramtools)

# gramtools
**TL;DR** Genome inference using prior information encoded as a refernce graph.

Gramtools builds a directed acyclical graph (DAG) of genetic variance from a population of genomes. Given a sample of reads, the graph is then annotated with quasimap coverage information. A personal reference genome for the population can then be infered from the annotated graph.

## Install
```pip3 install git+https://github.com/iqbal-lab-org/gramtools```

## Usage
Gramtools currently consists of three commands. These commands are documented in the wiki (see links below). In dependancy order, they are:
1) [build](https://github.com/iqbal-lab-org/gramtools/wiki/Commands%3A-build) - given a VCF and reference, construct the graph

2) [quasimap](https://github.com/iqbal-lab-org/gramtools/wiki/Commands%3A-quasimap) - given a set of reads and a graph, generate mapping coverage information

3) [infer](https://github.com/iqbal-lab-org/gramtools/wiki/Commands%3A-infer) - given coverage information and a grpah, infer a maximum likelihood genome

```
Gramtools

Usage:
  gramtools build --gram-directory --vcf --reference
  gramtools quasimap --gram-directory --reads --output-directory
  gramtools infer --gram-directory --mean-depth --error-rate --quasimap-directory --output
  gramtools (-h | --help)
  gramtools --version

Options:
  --gram-directory      Gramtools supporting data structures
  --quasimap-directory  Quasimap coverage information output directory
  --output-directory    Command output directory
  --output              Command output file path
  --mean-depth          Mean read coverage depth
  --error-rate          Per base error rate

  -h --help             Show this screen
  --version             Show version

Subcommands:
  build         Construct the graph and supporting data structures
  quasimap      Generate mapping coverage information given reads
  infer         Infer a maximum likelihood genome
```

Examples and documentation can be found in the GitHub [wiki](https://github.com/iqbal-lab-org/gramtools/wiki).


## License

MIT