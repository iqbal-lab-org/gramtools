import logging
import itertools
import collections
from Bio import SeqIO

from . import prg

log = logging.getLogger('gramtools')


def parse_args(common_parser, subparsers):
    parser = subparsers.add_parser('kmers',
                                   parents=[common_parser])

    parser.add_argument('--kmer-size', help='',
                        type=int,
                        required=True)
    parser.add_argument('--reference', help='',
                        type=str,
                        required=True)
    parser.add_argument('--output-fpath', help='',
                        type=str,
                        required=True)

    parser.add_argument('--kmer-region-size',
                        dest='kmer_region_size',
                        help='',
                        type=int)
    parser.add_argument('--all-kmers', help='',
                        action='store_true')


def _truncate_overshoot_alleles(max_overshoot_bases, inrange_alleles, reverse):
    if not inrange_alleles:
        return inrange_alleles

    if reverse:
        truncated_alleles = [allele[-max_overshoot_bases:]
                             for allele in inrange_alleles[0]]
        inrange_alleles[0] = truncated_alleles

    else:
        truncated_alleles = [allele[:max_overshoot_bases]
                             for allele in inrange_alleles[-1]]
        inrange_alleles[-1] = truncated_alleles
    return inrange_alleles


def _directional_alleles_range(kmer_region_size, start_region, regions, reverse):
    """Return list of regions which are within a given base distance
    of a starting regions.

    The starting region is not included in the returned region range.
    """
    min_distance = 0
    max_distance = 0
    last_max_allele_len = 0
    passed_kmer_region = False

    in_range_alleles = collections.deque()
    for region in regions.range(start_region, reverse=reverse):
        alleles = list(region.alleles)
        if reverse:
            in_range_alleles.appendleft(alleles)
        else:
            in_range_alleles.append(alleles)

        region_edge_distance = min_distance + region.min_allele_len
        passed_kmer_region = region_edge_distance >= kmer_region_size
        if passed_kmer_region:
            break

        min_distance += region.min_allele_len
        max_distance += region.max_allele_len
        last_max_allele_len = region.max_allele_len

    max_overshoot_bases = kmer_region_size - max_distance
    if not passed_kmer_region:
        max_overshoot_bases -= last_max_allele_len

    return in_range_alleles, max_overshoot_bases


def _alleles_within_range(kmer_region_size, start_region, regions):
    """Return regions within a max base distance (kmer_region_size)
    of start region. Start region is not included in the distance measure.
    """
    reverse_alleles, max_overshoot_bases = \
        _directional_alleles_range(kmer_region_size,
                                   start_region, regions,
                                   reverse=True)
    reverse_alleles = _truncate_overshoot_alleles(max_overshoot_bases,
                                                  reverse_alleles,
                                                  reverse=True)

    star_region_alleles = collections.deque([start_region.alleles])

    forward_alleles, max_overshoot_bases = \
        _directional_alleles_range(kmer_region_size,
                                   start_region, regions,
                                   reverse=False)
    forward_alleles = _truncate_overshoot_alleles(max_overshoot_bases,
                                                  forward_alleles,
                                                  reverse=False)

    return reverse_alleles + star_region_alleles + forward_alleles


def _genome_paths(regions):
    return itertools.product(*regions)


def _kmers_from_genome_path(path_parts, kmer_size):
    path = []
    for path_part in path_parts:
        path.extend(path_part)

    for i in range(len(path)):
        if i + kmer_size > len(path):
            break
        kmer = tuple(path[i: i + kmer_size])
        yield kmer


def _kmers_from_genome_paths(genome_paths, kmer_size, unique=False):
    """Generate all kmers which exist for a list of genome paths."""
    for path_parts in genome_paths:
        kmers = _kmers_from_genome_path(path_parts, kmer_size)

        seen_kmers = set()
        for kmer in kmers:
            if unique:
                if kmer not in seen_kmers:
                    seen_kmers.add(kmer)
                    yield kmer
            else:
                yield kmer


def _generate_variant_kmers(kmer_region_size, kmer_size, regions):
    """Generate kmers centered on variant sites."""
    seen_kmers = set()
    _variant_regions = (r for r in regions if r.is_variant_site)
    for start_region in _variant_regions:
        alleles = _alleles_within_range(kmer_region_size,
                                        start_region,
                                        regions)
        genome_paths = _genome_paths(alleles)
        kmers = _kmers_from_genome_paths(genome_paths, kmer_size,
                                         unique=True)
        for kmer in kmers:
            if kmer not in seen_kmers:
                yield kmer
                seen_kmers.add(kmer)


def _generate_all_kmers(kmer_size, regions):
    """Generate all kmers."""
    seen_kmers = set()
    for start_region in regions:
        kmer_region_size = kmer_size
        alleles, _ = _directional_alleles_range(kmer_region_size, start_region,
                                                regions, reverse=False)

        star_region_alleles = collections.deque([start_region.alleles])
        alleles = star_region_alleles + alleles

        genome_paths = _genome_paths(alleles)
        kmers = _kmers_from_genome_paths(genome_paths, kmer_size,
                                         unique=True)
        for kmer in kmers:
            if kmer not in seen_kmers:
                yield kmer
                seen_kmers.add(kmer)


def _generate(kmer_region_size, kmer_size, regions, all_kmers):
    """Generate kmers."""
    if all_kmers:
        generator = _generate_all_kmers(kmer_size, regions)
    else:
        generator = _generate_variant_kmers(kmer_region_size,
                                            kmer_size, regions)
    return (kmer for kmer in generator)


def _dump_kmers(kmers, output_fpath):
    count_written = 0
    with open(output_fpath, 'w') as fhandle:
        for kmer in kmers:
            fhandle.write(''.join(kmer) + '\n')
            count_written += 1
    return count_written


def _dump_masks(sites, alleles, args):
    log.debug('Writing sites mask to file')
    with open(args.sites_mask_fpath, 'w') as fhandle:
        fhandle.write(sites)
    log.debug('Writing alleles mask to file')
    with open(args.allele_mask_fpath, 'w') as fhandle:
        fhandle.write(alleles)


def run(args):
    log.info('Start process: generate kmers')

    log.debug('Parsing PRG')
    fasta_seq = str(SeqIO.read(args.reference, 'fasta').seq)
    regions = prg.parse(fasta_seq)

    log.debug('Generating Kmers')
    kmers = _generate(args.kmer_region_size,
                      args.kmer_size, regions,
                      args.all_kmers)

    log.debug('Writing kmers to file')
    count_written = _dump_kmers(kmers, args.output_fpath)
    log.debug('Number of kmers writen to file: %s', count_written)

    log.info('End process: generate kmers')
