# gramtools

Genome inference (and variant calling) with a reference graph of genetic variation, based on a vBWT.
A vBWT is a compressed suffix array used to encode both coordinates and known 
genetic variation.


# Build
You need C++11 and the boost library installed. Download with:

```
git clone --recursive https://github.com/iqbal-lab/gramtools
```

Run installation script:

```
./install.sh
```

Add the ``lib`` directory to your LD_LIBRARY_PATH using, e.g.

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/lib
```

Gramtools needs the boost/functional/hash.hpp header, so you may need to modify
the BOOST_HEADERS variable in the Makefile with the path to your boost
installation, something like path/to/boost_1_xx_x. Then you can compile:

```
make
```

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




This will create the following files which you do not need to look at:

1. out/csa - the compressed suffix array
2. out/bin - the PRG in integer encoding
3. out/memlog - a file used by SDSL for memory tracking
4. out/reads - a text file showing progress on how many reads have been processed.

and the following file which you will want to look at or parse:

1. out/covg - a text file with one line per site. Each line has a space-separated list of coverages (one per allele)


We'll be adding new scripts/code for processing.





