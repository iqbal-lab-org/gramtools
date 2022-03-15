[![Build Status](https://travis-ci.org/iqbal-lab-org/gramtools.svg?branch=master)](https://travis-ci.org/iqbal-lab-org/gramtools)
[![Docker Repository on Quay](https://quay.io/repository/iqballab/gramtools/status "Docker Repository on Quay")](https://quay.io/repository/iqballab/gramtools)

# gramtools
**TL;DR** genotype genetic variants using genome graphs.

`gramtools` builds a genome graph from a reference genome and a set of variants (VCF files or MSAs files). 
Given sequence data from an individual, the graph is annotated with coverage and genotyped, producing a VCF file, and a [jVCF file](https://github.com/iqbal-lab-org/jVCF-spec), of all the variation in the graph.
 
 A personalised reference genome for the sample is also inferred and new variation can be discovered 
 against it. You can then build a new graph from the initial and the new variants, and genotype this augmented graph.
 
 Check out [usage](#usage) for details, and [limitations](#limitations-and-recommendations) for 
 limitations and recommendations.

## Contents

- [Install](#installrun)
  - [Container](#from-container)
  - [Source](#from-source)
- [Usage](#usage)
- [Limitations](#limitations-and-recommendations)
- [Docs](#documentation)
- [Contribute](#contributing)
- [Licence](#licence)
- [Cite](#cite)

## Install/Run

### From container
The easiest way to run `gramtools` is via a container ([hosted on quay.io](https://quay.io/repository/iqballab/gramtools?tab=tags)).

To run with [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/):
```sh
tag="latest" # or, a specific released version
# Run with docker
docker run "quay.io/iqballab/gramtools:${tag}"
# Or run with singularity
URI="docker://quay.io/iqballab/gramtools:${tag}"
singularity exec "$URI" gramtools
```

### From source

#### Latest release

We recommend installing inside a virtual environment:
```sh
python -m venv venv_gramtools && source venv_gramtools/bin/activate
VERSION="1.10.0"
wget -O - "https://github.com/iqbal-lab-org/gramtools/releases/download/v${VERSION}/gramtools-${VERSION}.tar.gz" | tar xfz -
pip install ./gramtools-"${VERSION}"
```
gramtools contains a compiled backend component. The latest release includes a
precompiled backend for Linux.  Run `gramtools` to see if it is compatible with your
machine. If gramtools tells you it is not, run the following to compile the backend on
your machine:

```
pip install conan
bash ./gramtools-"${VERSION}"/build.sh
```

Please report any errors via the [issue tracker](https://github.com/iqbal-lab-org/gramtools/issues).

#### Latest source

This will always compile the binary.

```sh
git clone https://github.com/iqbal-lab-org/gramtools
pip install conan
bash ./gramtools/build.sh
pip install ./gramtools
```

#### Requirements

`Python >= 3.6` and `pip >= 20.0.2`.

If the backend needs to be compiled, you also need `CMake >= 3.1.2` and a C++17 
compatible compiler: `g++ >=8` (tested) or `clang >=7` (untested).

For `gramtools discover` to work, you additionally need to install `py-cortex-api` (it
is installable via PyPi/pip) and have `R` and `Perl` on your $PATH.

## Usage

```
$ gramtools -h

Usage: 
    gramtools [-h] [--debug] [--force] subcommand
    
    Subcommands:
        gramtools build -o GRAM_DIR --ref REFERENCE
                       (--vcf VCF [VCF ...] | --prgs_bed PRGS_BED | --prg PRG)
                       [--kmer_size KMER_SIZE] [--max_threads MAX_THREADS]

        gramtools genotype -i GRAM_DIR -o GENO_DIR
                          --reads READS [READS ...] --sample_id SAMPLE_ID
                          [--ploidy {haploid,diploid}]
                          [--max_threads MAX_THREADS] [--seed SEED]

        gramtools discover -i GENO_DIR -o DISCO_DIR
                          [--reads READS [READS ...]]

        gramtools simulate --prg PRG
                           [--max_num_paths MAX_NUM_PATHS]
                           [--sample_id SAMPLE_ID] [--output_dir OUTPUT_DIR]
```

### Subcommands explained
1)  [build][wiki_build] - 
    given a reference genome and VCF files or MSA files, builds an indexed graph. 

    For multiple-sequence alignment (MSA) files, graphs are 
    built by our tool [`make_prg`][make_prg]. You can provide the MSA files, or prg
    files you made using `make_prg`, using option `--prgs_bed`.
    To genotype complex regions (e.g. SNPs + SVs, or variants on multiple references),
    you must use this option.

    * `--prgs_bed`: builds a genome graphs of a set of regions
       with the rest of the linear reference genome in-between. The bed file 
       specifies each region to build (columns 1-3) and the name of a MSA or prg file correspoding 
       to the region in column 4. (see bed format [specification][bed])
       **Important**: currently the reference genome sequence must be the first entry in each MSA
       (also when using `make_prg` yourself)
    * `--kmer_size`: used for indexing the graph in preparation for
       `genotype`. higher `k` <=> faster `genotype`, but `build` output will consume more 
       disk space.
    * `--max_threads`: use this to build multiple regions in parallel when using
      `--prgs_bed`

2)  [genotype](https://github.com/iqbal-lab-org/gramtools/wiki/Commands%3A-genotype) - 
    map reads to a graph generated in `build` and genotype the graph. Produces genotype calls (VCF)
    and a personalised reference genome (fasta).
    * `--reads`: 1+ reads files in (fasta/fastq/sam/bam/cram) format
    * `--sample_id`: displayed in VCF & personalised reference outputs

3) [discover](https://github.com/iqbal-lab-org/gramtools/wiki/Commands%3A-discover) - 
discovers new variation against the personalised reference genome from `genotype` using
 one or more variant callers (currently: cortex).
 
4) simulate- samples paths through a prg, producing a fasta of the paths and a genotyped JSON
of the variant bubbles the path went through.
    * `--prg`: a prg file as output by `build`

## Limitations and recommendations

* gramtools is primarily a genotyper, not a variant caller. Variant discovery 
  with the `discover` command relies on existing tools
* gramtools is designed to use short, low-error rate reads only (e.g.
Illumina). We recommend trimming adapters off reads (e.g. using [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)) before genotyping with gramtools.
* gramtools currently performs exact matching of reads only, so relatively high read coverage (e.g.
  \>20X) is recommended
* gramtools does not currently genotype copy-number variants (CNVs) or samples with mixed ploidy

## Documentation

Examples, documentation, and planned future enhancements can be found in the [wiki](https://github.com/iqbal-lab-org/gramtools/wiki).

For the C++ source code, [doxygen](http://doxygen.nl/) formatted documentation can be generated by running 
```doxygen doc/Doxyfile.in```
from inside the gramtools directory.

The documentation gets generated in doc/html/index.html and provides a useful reference for all files, classes, functions and data structures in gramtools.

## Contributing

Contributions are always welcome!
Please refer to [CONTRIBUTING.md](CONTRIBUTING.md). For developer tips, head to the [wiki](https://github.com/iqbal-lab-org/gramtools/wiki/Developers%3A-tips).


## License

MIT

## Cite

If you use gramtools in your work, cite as:
> Letcher, B., Hunt, M. & Iqbal, Z. Gramtools enables multiscale variation analysis with genome graphs. *Genome Biology* 22, 259 (2021). https://doi.org/10.1186/s13059-021-02474-0


To cite the vBWT data structure: 
> Maciuca, S., del Ojo Elias, C., McVean, G., Iqbal, Z. A Natural Encoding of Genetic Variation in a Burrows-Wheeler Transform to Enable Mapping and Genome Inference. WABI2016: Algorithms in Bioinformatics. Lecture Notes in Computer Science, vol 9838. Springer, Cham. https://doi.org/10.1007/978-3-319-43681-4_18

To cite `make_prg` graph construction:
> Colquhoun, R.M., Hall, M.B., Lima, L. et al. Pandora: nucleotide-resolution bacterial pan-genomics with reference graphs. *Genome Biology* 22, 267 (2021). https://doi.org/10.1186/s13059-021-02473-1



[make_prg]: https://github.com/iqbal-lab-org/make_prg
[bed]: https://en.wikipedia.org/wiki/BED_(file_format)
[wiki_build]: https://github.com/iqbal-lab-org/gramtools/wiki/Commands%3A-build
