import logging
import itertools

from Bio import SeqIO

from . import genome_regions
from . import parse_prg


log = logging.getLogger('gramtools')


def _filter_regions(regions, nonvariant_kmers):
    """TODO"""
    if nonvariant_kmers:
        return iter(regions)

    for region in regions:
        if region.variant_site_marker is not None:
            yield region


def _directional_region_range(max_base_distance, start_region,
                              regions, reverse):
    """TODO"""
    tot_distance = 0
    regions_in_range = []
    for region in regions.range(start_region, reverse=reverse):
        if region == start_region:
            continue
        if tot_distance >= max_base_distance:
            break

        regions_in_range.append(region)
        tot_distance += region.min_allele_len

    if reverse:
        regions_in_range.reverse()
    return regions_in_range


def _regions_within_distance(max_base_distance, start_region,
                             regions):
    """TODO"""
    reverse_range = _directional_region_range(max_base_distance,
                                              start_region,
                                              regions,
                                              reverse=True)
    forward_range = _directional_region_range(max_base_distance,
                                              start_region,
                                              regions,
                                              reverse=False)
    return reverse_range + [start_region] + forward_range


def _genome_paths(region_range):
    """TODO"""
    range_alleles = [region.alleles for region in region_range]
    for genome_path in itertools.product(*range_alleles):
        yield ''.join(genome_path)


def _kmers_from_genome_paths(genome_paths, kmer_size):
    """TODO"""
    seen_kmers = set()
    for path in genome_paths:
        for i in range(len(path)):
            if i + kmer_size > len(path):
                break
            kmer = path[i: i + kmer_size]

            if kmer not in seen_kmers:
                seen_kmers.add(kmer)
                yield kmer


def generate(max_base_distance, kmer_size, regions, nonvariant_kmers):
    """TODO"""
    for start_region in _filter_regions(regions, nonvariant_kmers):
        region_range = _regions_within_distance(max_base_distance,
                                                start_region,
                                                regions)
        genome_paths = _genome_paths(region_range)
        kmers = _kmers_from_genome_paths(genome_paths, kmer_size)
        return kmers


def _dump_kmers(kmers, output_fpath):

    if output_fpath is None:
        log.debug('Writing kmers to stdout')
        for kmer in kmers:
            print(kmer)
        return

    log.debug('Writing kmers to file')
    with open(output_fpath, 'w') as fhandle:
        for kmer in kmers:
            fhandle.write(kmer + '\n')


def _dump_masks(sites, alleles):
    with open('mask_sites.txt', 'w') as fhandle:
        fhandle.write(sites)

    with open('mask_alleles.txt', 'w') as fhandle:
        fhandle.write(alleles)


def run(args):
    log.info('Start process: kmers')
    mask = True

    fasta_seq = str(SeqIO.read(args.reference, 'fasta').seq)
    regions = parse_prg.parse(fasta_seq)

    kmers = generate(args.kmer_region_distance,
                     args.kmer_size, regions,
                     args.nonvariant_kmers)
    _dump_kmers(kmers, args.output_fpath)

    if mask:
        sites_mask = genome_regions.sites_mask(regions)
        alleles_mask = genome_regions.alleles_mask(regions)
        _dump_masks(sites_mask, alleles_mask)
    log.info('End process: build')
