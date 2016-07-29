# gramtools

Mapping and variant calling with a reference graph of genetic variation, based on a vBWT.
A vBWT is a compressed suffix array used to encode both coordinates and known 
genetic variation.

gramtools is very much a work in progress. We've so far proven that we can encode
genetic variation in a BWT in such a way as to enable variation-aware mapping, and so
that the notion of allelism is preserved. When there is ambiguity in mapping, the structure
naturally reveals if the alternate mapping options are paralogs or alleles.

We've implemented backward-search, the infrastructure for bidirectional search (so we can find MEMS),
exact matching. For a haploid organism, this enables us to find a reference genome which is a mosaic
of the sequence in the graph, naturally recapitulating what we expect for biological reasons. Having found 
this nearby reference, traditional mappign and variant calling will work much better.

In terms of RAM usage gramtools scales to human genomes with millions of SNPs/indels in the graph, but
currently mapping speed is too slow on human-size graphs. We have a list of very promising performance
upgrades coming soon though...

Forthcoming changes:
 - performance improvements (eg through storing ranks)
 - inexact matching
 - phasing
 - inference of a pair/tuple of haplotypes using an HMM as in Dilthey et al (2015)
 - full MEM-based alignment to the inferred nearby haplotypes.

Do read our paper, which you can find here: 

NB we have an open bug at the moment where, after the program has finished, right at the very end, 
when the program finishes the shutdown process causes a crash. We're looking into it!

Sorina, Carlos, Gil, Zam

# Compiling
You need C++11. Go into the src/ directory and type `make`.

# Usage

* Build the PRG file from a VCF
```
perl utils/vcf_to_linear_prg.pl --outfile my_species.prg --vcf big.vcf --ref reference.fa
```
This will create my_species.prg, my_species.prg.fa, my_species.prg.mask_alleles and my_species.prg.mask_sites

* Create a kmers file

```
python utils/variantKmers.py -f my_species.prg.fa  -k 9 -n > my_species.kmers.txt
```

* Map reads to the PRG, putting output in directory out/

```
gramtools --prg my_species.prg --csa out/csa --ps my_species.prg.mask_sites --pa my_species.prg.mask_alleles
          --co out/covg --ro out/reads --po out/bin --log out/memlog --ksize 9 --kfile my_species.kmers.txt 
          --input sample.fastq
```






This will create
