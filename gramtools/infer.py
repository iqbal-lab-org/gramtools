import json
import logging

import vcf

from . import genotyper
from . import paths

log = logging.getLogger('gramtools')
PRG_READ_CHUNK_SIZE = 10000
DIGITS = {
    '0', '1', '2', '3',
    '4', '5', '6', '7',
    '8', '9'
}


def parse_args(common_parser, subparsers):
    parser = subparsers.add_parser('infer',
                                   parents=[common_parser])
    parser.add_argument('--gram-directory',
                        help='',
                        type=str,
                        required=True)
    parser.add_argument('--quasimap-directory',
                        help='',
                        type=str,
                        required=True)
    parser.add_argument('--mean-depth',
                        help='',
                        type=int,
                        required=True)
    parser.add_argument('--error-rate',
                        help='',
                        type=float,
                        required=True)
    parser.add_argument('--haploid',
                        help='',
                        action='store_true',
                        required=False)
    parser.add_argument('--output-fasta',
                        help='',
                        type=str,
                        required=False)
    parser.add_argument('--output-vcf',
                        help='',
                        type=str,
                        required=False)


def run(args):
    log.info('Start process: infer')
    _paths = paths.generate_infer_paths(args)

    all_per_base_coverage = _load_per_base_coverage(_paths['allele_base_coverage'])
    allele_groups, all_groups_site_counts = \
        _load_grouped_allele_coverage(_paths['grouped_allele_counts_coverage'])

    allele_indexes = _max_likelihood_alleles(args.mean_depth,
                                             args.error_rate,
                                             all_per_base_coverage,
                                             all_groups_site_counts,
                                             allele_groups)

    prg_parser = _parse_prg(_paths['prg'])

    if args.output_fasta is not None:
        writer = _FastaWriter(args.output_fasta)
        _dump_fasta(prg_parser, allele_indexes, writer)
        writer.close()
    elif args.output_vcf is not None:
        _dump_vcf(allele_indexes, _paths, args)
    else:
        log.error("Either --output-fasta or --output-vcf must be specified")
        exit(1)

    log.info('End process: infer')


def _dump_vcf(allele_indexes, _paths, args):
    if _paths['perl_generated_vcf'] is '' or _paths['perl_generated_vcf'] is None:
        log.error('VCF file must be used as part of gramtools build command. Exiting.')
    
    file_handle = open(_paths['perl_generated_vcf'], 'r')
    vcf_reader = vcf.Reader(file_handle)

    record_attributes = [
        'CHROM',
        'POS',
        'ID',
        'REF',
        'ALT',
        'QUAL',
        'FILTER',
        'INFO',
        'FORMAT',
        '_sample_indexes'
    ]

    with open(args.output_vcf, 'w') as file_handle:
        vcf_writer = vcf.Writer(file_handle, vcf_reader)

        for record, allele_index in zip(vcf_reader, allele_indexes):
            if allele_index is None or allele_index == 0:
                continue

            attributes = {name: getattr(record, name) for name in record_attributes}
            attributes['ALT'] = [attributes['ALT'][allele_index - 1]]
            attributes = [attributes[name] for name in record_attributes]
            new_record = vcf.model._Record(*attributes)

            vcf_writer.write_record(new_record)


def _dump_fasta(prg_parser, allele_indexes, writer):
    allele_index = next(allele_indexes)

    for cursor in prg_parser:
        if cursor.just_left_site:
            try:
                allele_index = next(allele_indexes)
            except StopIteration:
                allele_index = None

            # no coverage data, use first allele
            if allele_index is None:
                allele_index = 0

        if cursor.on_marker:
            continue

        if not cursor.within_allele:
            writer.append(cursor.char)
            continue

        if cursor.allele_id == allele_index:
            writer.append(cursor.char)


class _FastaWriter:
    def __init__(self, fpath):
        self._fpath = fpath
        self._fhandle = open(self._fpath, 'w')
        self._fhandle.write('> A personal reference sequence built using gramtools.\n')

        self._cache = []
        self._max_cache_size = 10000

    def append(self, char):
        self._cache.append(char)
        if len(self._cache) > self._max_cache_size:
            self._flush()

    def _flush(self):
        cache = ''.join(self._cache)
        self._fhandle.write(cache)
        self._cache = []

    def close(self):
        self._flush()
        self._fhandle.close()


def _max_likelihood_alleles(mean_depth,
                            error_rate,
                            all_per_base_coverage,
                            all_groups_site_counts,
                            allele_groups):
    sites_coverage = iter(zip(all_per_base_coverage, all_groups_site_counts))
    for per_base_coverage, groups_site_counts in sites_coverage:
        allele_index = _max_likelihood_allele(mean_depth,
                                              error_rate,
                                              per_base_coverage,
                                              groups_site_counts,
                                              allele_groups)
        yield allele_index


def _max_likelihood_allele(mean_depth,
                           error_rate,
                           per_base_coverage,
                           groups_site_counts,
                           allele_groups):
    gtyper = genotyper.Genotyper(mean_depth,
                                 error_rate,
                                 groups_site_counts,
                                 per_base_coverage,
                                 allele_groups)
    gtyper.run()

    if gtyper.likelihoods is None:
        return None

    for alleles, likelihood in gtyper.likelihoods:
        if len(alleles) != 1:
            continue
        return alleles.pop()
    raise ValueError('Genotyper likelihoods does not have any valid haploid data')


def _load_grouped_allele_coverage(fpath):
    with open(fpath, 'r') as fhandle:
        data = json.load(fhandle)
    groups_coverage = data['grouped_allele_counts']

    allele_groups = groups_coverage['allele_groups']
    for key, value in allele_groups.items():
        allele_groups[key] = set(value)

    groups_site_counts = groups_coverage['site_counts']
    return allele_groups, groups_site_counts


def _load_per_base_coverage(fpath):
    with open(fpath, 'r') as fhandle:
        data = json.load(fhandle)
    data = data['allele_base_counts']
    return data


def _is_int(data):
    if data is None or data == '':
        return False
    if len(data) > 1:
        return True
    return data in DIGITS


def _parse_prg_chars(chars):
    int_chars = []
    for char in chars:
        if _is_int(char):
            int_chars.append(char)
            continue

        else:
            if int_chars:
                yield ''.join(int_chars)
                int_chars = []
            yield char

    if int_chars:
        yield ''.join(int_chars)


class _Cursor:
    char = None
    site_marker = None
    allele_id_counter = None

    on_marker = False
    just_left_site = False

    @property
    def allele_id(self):
        if not self.within_allele:
            return None
        return self.allele_id_counter

    @property
    def within_allele(self):
        return self.site_marker is not None and self.on_marker is False


def _parse_prg_structure(chars, cursor):
    for char in _parse_prg_chars(chars):
        cursor.char = char

        if cursor.just_left_site:
            cursor.just_left_site = False
            cursor.site_marker = None
            cursor.allele_id_counter = None

        if not _is_int(char):
            cursor.on_marker = False
            yield cursor
            continue

        cursor.on_marker = True

        site_boundary = int(char) % 2 != 0
        if site_boundary:
            entring = cursor.site_marker is None
            if entring:
                cursor.site_marker = int(char)
                cursor.allele_id_counter = 0
                yield cursor
                continue

            exiting = cursor.site_marker is not None
            if exiting:
                cursor.just_left_site = True
                cursor.allele_id_counter = None
                yield cursor
                continue

        allele_boundary = int(char) % 2 == 0
        if allele_boundary:
            cursor.allele_id_counter += 1
            yield cursor


def _read_chunk(file_handle, chunk_size=PRG_READ_CHUNK_SIZE):
    chars = file_handle.read(chunk_size)
    if not chars:
        return None

    if _is_int(chars[-1]):
        extra_chars = []
        while not extra_chars or _is_int(extra_chars[-1]):
            char = file_handle.read(1)
            if not chars:
                break
            extra_chars.append(char)
        chars += ''.join(extra_chars)
    return chars


def _parse_prg(prg_fpath):
    cursor = _Cursor()
    with open(prg_fpath) as file_handle:
        while True:
            chars = _read_chunk(file_handle)
            if chars is None:
                break

            for cursor in _parse_prg_structure(chars, cursor):
                yield cursor
