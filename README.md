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
