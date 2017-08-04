import logging
import itertools

from Bio import SeqIO

from . import genome_regions
from . import parse_prg

log = logging.getLogger('gramtools')


def parse_args(common_parser, subparsers):
    parser = subparsers.add_parser('kmers',
                                   parents=[common_parser])
    parser.add_argument('--kmer-size', help='',
                        type=int)
    parser.add_argument('--reference', help='',
                        type=str)
    parser.add_argument('--kmer-region-distance',
                        dest='kmer_region_distance',
                        help='',
                        type=int)
    parser.add_argument('--nonvariant-kmers', help='',
                        action='store_true')
    parser.add_argument('--output-fpath', help='',
                        type=str)


def _filter_regions(regions, nonvariant_kmers):
    """Yield regions with filtering based on a region's variant site status."""
    if nonvariant_kmers:
        for region in regions:
            yield region
        return

    for region in regions:
        if region.is_variant_site:
            yield region


def _directional_region_range(max_base_distance, start_region,
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
        if tot_distance >= max_base_distance:
            break

        regions_in_range.append(region)
        tot_distance += region.min_allele_len

    if reverse:
        regions_in_range.reverse()
    return regions_in_range


def _regions_within_distance(max_base_distance, start_region,
                             regions):
    """Return regions within a max base distance of start region.
    Start region is not included in the distance measure.
    """
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


def _generate(max_base_distance, kmer_size, regions, nonvariant_kmers):
    """Generate kmers."""
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

    else:
        log.debug('Writing kmers to file')
        with open(output_fpath, 'w') as fhandle:
            for kmer in kmers:
                fhandle.write(kmer + '\n')


def _dump_masks(sites, alleles):
    log.debug('Writing sites and alleles masks to file')
    with open('mask_sites.txt', 'w') as fhandle:
        fhandle.write(sites)
    with open('mask_alleles.txt', 'w') as fhandle:
        fhandle.write(alleles)


def run(args):
    log.info('Start process: generate kmers')
    mask = True

    fasta_seq = str(SeqIO.read(args.reference, 'fasta').seq)
    regions = parse_prg.parse(fasta_seq)

    kmers = _generate(args.kmer_region_distance,
                      args.kmer_size, regions,
                      args.nonvariant_kmers)
    _dump_kmers(kmers, args.output_fpath)

    if mask:
        sites_mask = genome_regions.sites_mask(regions)
        alleles_mask = genome_regions.alleles_mask(regions)
        _dump_masks(sites_mask, alleles_mask)
    log.info('End process: generate kmers')
