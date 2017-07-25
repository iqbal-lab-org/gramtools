import argparse

from Bio import SeqIO

from . import genome_regions
from . import generate_kmers
from . import parse_prg
from . import random_genome


def _dump_kmers(kmers, output_fpath):
    if output_fpath is None:
        for kmer in kmers:
            print(kmer)
        return

    with open(output_fpath, 'w') as fhandle:
        for kmer in kmers:
            fhandle.write(kmer + '\n')


def _dump_masks(sites, alleles):
    with open('mask_sites.txt', 'w') as fhandle:
        fhandle.write(sites)

    with open('mask_alleles.txt', 'w') as fhandle:
        fhandle.write(alleles)


def run(fasta_fpath, kmer_size, nonvariant_kmers, num_reads,
        mask, output_fpath=None):
    fasta_seq = str(SeqIO.read(fasta_fpath, 'fasta').seq)
    regions = parse_prg.parse(fasta_seq)

    max_base_distance = 20
    kmers = generate_kmers.generate_kmers(max_base_distance,
                                          kmer_size, regions,
                                          nonvariant_kmers)
    _dump_kmers(kmers, output_fpath)

    if mask:
        sites_mask = genome_regions.sites_mask(regions)
        alleles_mask = genome_regions.alleles_mask(regions)
        _dump_masks(sites_mask, alleles_mask)

    if num_reads:
        random_genome.dump_random_genome(num_reads, regions)


if __name__ == "__main__":
    aparser = argparse.ArgumentParser(description='PRG kmer generator')
    aparser.add_argument("-f", dest="fasta", help="Fasta file", required=True)
    aparser.add_argument("-n", dest="nonvariant_kmers", default=False,
                         help="Print kmers from intervariant _regions",
                         action='store_true')
    aparser.add_argument("-k", dest="ksize", help="kmer size [31]",
                         type=int, default=31)
    aparser.add_argument("-r", dest="nreads",
                         help="Generates random genome and random reads",
                         type=int, default=0)
    aparser.add_argument("-m", dest="mask", action='store_true',
                         default=False, help="Dump mask (mask.txt)")
    options = aparser.parse_args()

    run(options.fasta, options.ksize, options.nonvariant_kmers,
        options.nreads, options.mask)
