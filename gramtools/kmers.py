import logging
import itertools

from Bio import SeqIO

from . import genome_regions
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
    parser.add_argument('--kmer-region-size',
                        dest='kmer_region_size',
                        help='',
                        type=int,
                        required=True)
    parser.add_argument('--output-fpath', help='',
                        type=str,
                        required=True)

    parser.add_argument('--nonvariant-kmers', help='',
                        action='store_true')


def _filter_regions(regions, nonvariant_kmers):
    """Yield regions with filtering based on a region's variant site status."""
    if nonvariant_kmers:
        for region in regions:
            yield region
        return

    for region in regions:
        if region.is_variant_site:
            yield region


def _directional_region_range(kmer_region_size, start_region,
                              regions, reverse):
    """Return list of regions which are within a given base distance
    of a starting regions.

    The starting region is not included in the returned region range.
    """
    tot_distance = 0
    regions_in_range = []
    for region in regions.range(start_region, reverse=reverse):
        if region == start_region:
            continue
        if tot_distance >= kmer_region_size:
            break

        regions_in_range.append(region)
        tot_distance += region.min_allele_len

    if reverse:
        regions_in_range.reverse()
    return regions_in_range


def _regions_within_distance(kmer_region_size, start_region,
                             regions):
    """Return regions within a max base distance of start region.
    Start region is not included in the distance measure.
    """
    reverse_range = _directional_region_range(kmer_region_size,
                                              start_region,
                                              regions,
                                              reverse=True)
    forward_range = _directional_region_range(kmer_region_size,
                                              start_region,
                                              regions,
                                              reverse=False)
    return reverse_range + [start_region] + forward_range


def _genome_paths(region_range):
    """Generate all genome paths which exist for a given region range."""
    range_alleles = [region.alleles for region in region_range]
    for genome_path in itertools.product(*range_alleles):
        yield ''.join(genome_path)


def _kmers_from_genome_paths(genome_paths, kmer_size):
    """Generate all kmers which exist for a list of genome paths."""
    seen_kmers = set()
    for path in genome_paths:
        for i in range(len(path)):
            if i + kmer_size > len(path):
                break
            kmer = path[i: i + kmer_size]

            if kmer not in seen_kmers:
                seen_kmers.add(kmer)
                yield kmer


def _generate(kmer_region_size, kmer_size, regions, nonvariant_kmers):
    """Generate kmers."""
    seen_kmers = set()
    for start_region in _filter_regions(regions, nonvariant_kmers):
        region_range = _regions_within_distance(kmer_region_size,
                                                start_region,
                                                regions)
        genome_paths = _genome_paths(region_range)
        kmers = _kmers_from_genome_paths(genome_paths, kmer_size)
        for kmer in kmers:
            if kmer not in seen_kmers:
                yield kmer
                seen_kmers.add(kmer)


def _dump_kmers(kmers, output_fpath):
    count_written = 0
    with open(output_fpath, 'w') as fhandle:
        for kmer in kmers:
            fhandle.write(kmer + '\n')
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
                      args.nonvariant_kmers)

    log.debug('Writing kmers to file')
    count_written = _dump_kmers(kmers, args.output_fpath)
    log.debug('Number of kmers writen to file: %s', count_written)

    if hasattr(args, 'sites_mask_fpath') and hasattr(args, 'allele_mask_fpath'):
        sites_mask = genome_regions.sites_mask(regions)
        alleles_mask = genome_regions.alleles_mask(regions)
        _dump_masks(sites_mask, alleles_mask, args)

    log.info('End process: generate kmers')
