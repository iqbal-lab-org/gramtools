import logging
import itertools
import collections

from . import prg_regions_parser

log = logging.getLogger('gramtools')


def setup_command_parser(common_parser, subparsers):
    parser = subparsers.add_parser('simulate',
                                   parents=[common_parser])
    parser.add_argument('--max-num-reads', help='',
                        type=int, default=None)
    parser.add_argument('--read-length', help='',
                        type=int)
    parser.add_argument('--reference', help='',
                        type=str)
    parser.add_argument('--output-fpath', help='',
                        type=str)


def _read_regions(read_length, start_region, regions):
    read_regions = collections.deque()

    bases_count = start_region.min_allele_len
    for region in regions.range(start_region, reverse=True):
        read_regions.appendleft(region)
        if bases_count >= read_length:
            break
        bases_count += region.min_allele_len

    read_regions.append(start_region)
    return read_regions


def _variants_read_regions(read_length, genome_regions):
    variant_regions = (region for region in genome_regions
                       if region.is_variant_site)

    for start_region in variant_regions:
        read_regions = _read_regions(read_length, start_region,
                                     genome_regions)
        yield read_regions


def _generate_genome_paths(read_regions):
    for regions in read_regions:
        all_alleles = (region.alleles for region in regions)
        for alleles in itertools.product(*all_alleles):
            genome_path = []
            for allele in alleles:
                for base in allele:
                    genome_path.append(base)
            yield tuple(genome_path)


def _generate_reads(read_length, genome_regions, max_num_reads=None):
    reads = set()
    all_read_regions = _variants_read_regions(read_length, genome_regions)
    for genome_path in _generate_genome_paths(all_read_regions):
        if max_num_reads is not None and len(reads) == max_num_reads:
            break

        read = genome_path[-read_length:]
        if len(read) != read_length:
            continue
        if read not in reads:
            reads.add(read)
            yield ''.join(read)
    log.debug('Number of reads generated: %s', len(reads))


def _dump_random_reads(reads, quality, output_fpath):
    with open(output_fpath, 'w') as fhandle:
        for i, read in enumerate(reads):
            fastq_line = '@{read_num}\n{read}\n+\n{quality}\n'.format(
                read_num=i, read=read, quality=quality)
            fhandle.write(fastq_line)


def run(args):
    log.info('Start process: simulate')

    log.debug('Parsing PRG')
    with open(args.reference, 'r') as file_handle:
        prg_seq = file_handle.read()

    regions = prg_regions_parser.parse(prg_seq)

    log.debug('Generating reads')
    reads = _generate_reads(args.read_length, regions, args.max_num_reads)
    read_qualities = 'H' * args.read_length
    _dump_random_reads(reads, read_qualities, args.output_fpath)

    log.info('End process: simulate')
