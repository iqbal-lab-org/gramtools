[![Build Status](https://travis-ci.org/iqbal-lab-org/gramtools.svg?branch=dev)](https://travis-ci.org/iqbal-lab-org/gramtools)

# gramtools
Genomes evolve by recombination and mutation. 
gramtools uses a reference graph to  model new genomes as recombinants (mosaics) of previous genomes,
(by "quasimapping" reads to the graph) and then genotyping within the graph.

## Motivation
gramtools finds the nearest mosaic from a reference panel of genomes to a sample(a "personalised reference"). The user can then either
* use traditional mapping-based variant calling, 
* or, wait until we finish implementing variant discovery on top of the graph.

We recommend the first option for now.

## Install
With root:

```sudo pip3 install git+https://github.com/iqbal-lab-org/gramtools```

Without root:

```python3 -m venv gramtools_virtualenv && source ./gramtools_virtualenv/bin/activate && pip3 install git+https://github.com/iqbal-lab-org/gramtools```

## Run
Initial step done just once (per species)
* Generate a graph from known genetic variation and a standard reference genome (gramtools does this automatically from a VCF, but you can generate your own from a multiple sequence alignment or some combination of both).

For each sample
* Quality trim the reads 
* Quasimap reads to the graph, and infer mosaic reference
* Use whatever standard tools you like to call variants (samtools, GATK, etc)
* Convert variants back to standard coordinates (not supported yet)

### Build Graph
```gramtools build --gram-directory ./gram --vcf ./vcf --reference ./reference --max-read-length 150```

| parameter           | description                                                     |
|---------------------|-----------------------------------------------------------------|
| `--gram-directory`  | output directory for gramtools build files (created if missing) |
| `--vcf`             | input variant call file detailing genetic variants              |
| `--reference`       | generic reference genome which compliments the VCF              |
| `--max-read-length` | maximum read length used during the `quasimap` command          |

### Quasimap reads to graph
```gramtools quasimap --gram-directory ./gram --reads ./reads```

| parameter          | description                                                              |
|--------------------|--------------------------------------------------------------------------|
| `--gram-directory` | output directory for gramtools build files (created if missing)          |
| `--reads`          | input read samples fastq file, this parameter can be used multiple times |
