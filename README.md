[![Build Status](https://travis-ci.org/iqbal-lab-org/gramtools.svg?branch=dev)](https://travis-ci.org/iqbal-lab-org/gramtools)

# gramtools
Infer a reference genome from a population of genomes for better read mapping coverage.

## Install
```sudo pip install git+https://github.com/iqbal-lab-org/gramtools```

## Running
Infering a reference genome from a population takes two steps:
1) Generate a graph from population genetic variants (VCF) and a generic reference.
2) Infer a reference genome from the graph by analysing read coverage with exact matching.

### Build graph
```gramtools build --gram-directory ./gram --vcf ./vcf --reference ./reference --max-read-length 150```

| parameter           | description                                                     |
|---------------------|-----------------------------------------------------------------|
| `--gram-directory`  | output directory for gramtools build files (created if missing) |
| `--vcf`             | input variant call file detailing genetic variants              |
| `--reference`       | generic reference genome which compliments the VCF              |
| `--max-read-length` | maximum read length used during the `quasimap` command          |

### Infer reference genome
```gramtools quasimap --gram-directory ./gram --reads ./reads```

| parameter          | description                                                     |
|--------------------|-----------------------------------------------------------------|
| `--gram-directory` | output directory for gramtools build files (created if missing) |
| `--reads`          | input read samples fastq file                                   |
