import random
import logging
import itertools
import collections

from . import prg

log = logging.getLogger('gramtools')


def parse_args(common_parser, subparsers):
    parser = subparsers.add_parser('simulate',
                                   parents=[common_parser])
    parser.add_argument('--max-num-reads', help='',
                        type=int, default=None)
    parser.add_argument('--read-length', help='', type=int)
    parser.add_argument('--prg', help='', type=str)
    parser.add_argument('--output-fpath', help='', type=str)
    parser.add_argument('--noise', help='', type=float, required=False)
    parser.add_argument('--variant-sites-only',
                        help='',
                        action='store_true',
                        default=False,
                        required=False)


def _read_regions(start_region, read_length, regions):
    read_regions = collections.deque()
    bases_count = start_region.min_allele_len
    for region in regions.range(start_region, reverse=True):
        read_regions.appendleft(region)
        if bases_count >= read_length:
            break
        bases_count += region.min_allele_len
    read_regions.append(start_region)
    return read_regions


def _generate_genome_paths(read_regions):
    for regions in read_regions:
        all_alleles = (r.alleles for r in regions)
        for alleles in itertools.product(*all_alleles):
            genome_path = []
            for allele in alleles:
                genome_path += list(allele)
            yield tuple(genome_path)


def _variant_read_regions(read_length, regions):
    variant_regions = (r for r in regions if r.is_variant_site)
    for start_region in variant_regions:
        read_regions = _read_regions(start_region, read_length, regions)
        yield read_regions


def _generate_variant_site_reads(read_length, regions, max_num_reads):
    previous_reads = set()
    read_regions = _variant_read_regions(read_length, regions)
    for genome_path in _generate_genome_paths(read_regions):
        if len(previous_reads) == max_num_reads:
            break

        read = genome_path[-read_length:]
        if len(read) != read_length:
            continue
        if read in previous_reads:
            continue

        previous_reads.add(read)
        yield ''.join(read)


def _random_read_regions(read_length, regions):
    for start_region in regions.iter_random():
        read_regions = _read_regions(start_region, read_length, regions)
        yield read_regions


def _generate_random_reads(read_length, regions, max_num_reads):
    previous_reads = set()
    read_regions = _random_read_regions(read_length, regions)
    for genome_path in _generate_genome_paths(read_regions):
        if len(previous_reads) == max_num_reads:
            break

        upper_read_start_index = len(genome_path) - read_length
        read_start_index = random.randint(0, upper_read_start_index)
        read = genome_path[read_start_index:]

        if read in previous_reads:
            continue
        previous_reads.add(read)
        yield ''.join(read)


def _dump_reads(reads, output_fpath):
    with open(output_fpath, 'w') as fhandle:
        for i, read in enumerate(reads):
            quality = 'H' * len(read)
            fastq_line = '@{read_num}\n{read}\n+\n{quality}\n'.format(
                read_num=i, read=read, quality=quality)
            fhandle.write(fastq_line)


def _add_noise(reads, noise_rate):
    other_bases = {
        'A': ['C', 'G', 'T'],
        'C': ['A', 'G', 'T'],
        'G': ['A', 'C', 'T'],
        'T': ['A', 'C', 'G'],
    }
    for read in reads:
        read = list(read)
        for i, base in enumerate(read):
            implement_noise = random.random() < noise_rate
            if not implement_noise:
                continue
            choice_bases = other_bases[base]
            read[i] = random.choice(choice_bases)
        read = ''.join(read)
        yield read


def run(args):
    log.info('Start process: simulate')

    log.debug('Parsing PRG')
    with open(args.prg, 'r') as file_handle:
        prg_seq = file_handle.read()

    regions = prg.parse(prg_seq)

    if args.variant_sites_only:
        log.debug('Generating reads overlapping variant sites only')
        reads = _generate_variant_site_reads(args.read_length,
                                             regions,
                                             args.max_num_reads)
    else:
        log.debug('Generating reads')
        reads = _generate_random_reads(args.read_length,
                                       regions,
                                       args.max_num_reads)

    if hasattr(args, 'noise'):
        log.debug('Adding noise to reads')
        reads = _add_noise(reads, args.noise)

    _dump_reads(reads, args.output_fpath)
    log.info('End process: simulate')
