import os
import shutil
import typing
import bisect
import logging

import vcf
import cortex
from Bio import SeqIO, Seq
import cluster_vcf_records

from .. import paths
from .. import prg

log = logging.getLogger('gramtools')


def parse_args(common_parser, subparsers):
    parser = subparsers.add_parser('discover',
                                   parents=[common_parser])
    parser.add_argument('--gram-directory',
                        help='',
                        type=str,
                        required=True)
    parser.add_argument('--inferred-vcf',
                        help='',
                        type=str,
                        required=True)
    parser.add_argument('--reads',
                        help='',
                        type=str,
                        required=True)
    parser.add_argument('--output-vcf',
                        help='',
                        type=str,
                        required=False)


def run(args):
    log.info('Start process: discover')
    _paths = paths.generate_discover_paths(args)
    os.mkdir(_paths['tmp_directory'])

    with open(_paths['prg'], 'r') as file_handle:
        prg_seq = file_handle.read()
    reference = _get_reference(prg_seq)
    fasta_description = "PRG base reference for gramtools"
    _dump_fasta(reference, fasta_description, _paths['base_reference'])

    file_handle = open(args.inferred_vcf, 'r')
    vcf_reader = vcf.Reader(file_handle)

    inferred_reference = _get_inferred_reference(reference, vcf_reader)
    fasta_description = "inferred personal reference generated using gramtools"
    _dump_fasta(inferred_reference, fasta_description, _paths['inferred_reference'])

    cortex.calls(_paths['inferred_reference'], args.reads, _paths['cortex_vcf'])

    inferred_reference_length = len(inferred_reference)
    rebased_vcf_records = _rebase_vcf(args.inferred_vcf,
                                      inferred_reference_length,
                                      _paths['cortex_vcf'])

    template_vcf_file_path = _paths['cortex_vcf']
    _dump_rebased_vcf(rebased_vcf_records, _paths['rebase_vcf'], template_vcf_file_path)

    vcf_files = [_paths['rebase_vcf']]
    cluster = cluster_vcf_records.vcf_clusterer.VcfClusterer(vcf_files, _paths['base_reference'], args.output_vcf)
    cluster.run()

    shutil.rmtree(_paths['tmp_directory'])
    log.info('End process: discover')


def _rebase_vcf(base_vcf_file_path,
                secondary_reference_length: int,
                secondary_vcf_file_path):
    """Rebase secondary_vcf so that it uses same reference as base_vcf.

    base_vcf -> VCF describing inferred haploid
    secondary_vcf -> VCF describing new variants with inferred haploid as reference
    """

    base_records = _load_records(base_vcf_file_path)
    secondary_records = _load_records(secondary_vcf_file_path)

    # Given a site index position in secondary_vcf, need to know if within base_vcf site
    secondary_regions = _get_secondary_regions(base_records, secondary_reference_length)

    new_vcf_records = [_rebase_vcf_record(vcf_record, secondary_regions)
                       for vcf_record in secondary_records]
    return new_vcf_records


def _dump_rebased_vcf(records, rebase_file_path, template_vcf_file_path):
    template_handle = open(template_vcf_file_path, 'r')
    template_vcf = vcf.Reader(template_handle)

    rebase_handle = open(rebase_file_path, 'w')
    writer = vcf.Writer(rebase_handle, template_vcf)
    for record in records:
        writer.write_record(record)
    writer.close()


class _Region:
    def __init__(self, start, end, offset, record):
        self.start = start
        self.end = end
        self.offset = offset
        self.record = record

    @property
    def is_site(self):
        return self.record is not None

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __repr__(self):
        return str(self.__dict__)

    def __lt__(self, other):
        if isinstance(other, _Region):
            return self.start < other.start
        else:
            return self.start < other

    def __len__(self):
        return self.end - self.start + 1


_Regions = typing.List[_Region]


def _get_secondary_regions(base_records, secondary_reference_length) -> _Regions:
    secondary_site_regions = _secondary_regions_for_base_sites(base_records)
    secondary_nonsite_regions = _secondary_regions_for_base_nonsites(secondary_site_regions,
                                                                     base_records,
                                                                     secondary_reference_length)

    secondary_regions = _merge_regions(secondary_site_regions, secondary_nonsite_regions)
    return secondary_regions


_vcf_record_attributes = [
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


def _make_vcf_record(**attributes):
    attributes['ALT'] = [vcf.model._Substitution(x) for x in attributes['ALT']]
    attributes = [attributes.get(name) for name in _vcf_record_attributes]

    new_vcf_record = vcf.model._Record(*attributes)
    return new_vcf_record


def _modify_vcf_record(vcf_record, **new_attributes):
    current_attributes = {name: getattr(vcf_record, name, None)
                          for name in _vcf_record_attributes}
    attributes = current_attributes
    attributes.update(new_attributes)

    new_vcf_record = _make_vcf_record(**attributes)
    return new_vcf_record


def _rebase_vcf_record(vcf_record, secondary_regions: _Regions):
    first_region_index = _find_start_region_index(vcf_record, secondary_regions)
    overlap_regions = _find_overlap_regions(vcf_record, first_region_index, secondary_regions)

    # TODO: what about start and end in nonsites and sites in the middle?
    # below condition probably not general enough
    overlap_sites = any(region.is_site for region in overlap_regions)
    if not overlap_sites:
        offset_sum = sum(region.offset for region in overlap_regions)
        new_pos = vcf_record.POS + offset_sum

        vcf_record = _modify_vcf_record(vcf_record, POS=new_pos)
        return vcf_record

    if overlap_regions[0].is_site:
        vcf_record = _pad_vcf_record_start(vcf_record, overlap_regions)

    if overlap_regions[-1].is_site:
        vcf_record = _pad_vcf_record_end(vcf_record, overlap_regions)

    return vcf_record


def _pad_vcf_record_start(vcf_record, overlap_regions: _Regions):
    """Pad start of VCF record.

    Set REF to base site REF for given region. Update POS to region start.
    Pad ALT start with region start ALT (if any from current start backwards)
    """
    first_region = overlap_regions[0]

    record_start_index = vcf_record.POS - 1
    region_start_index = first_region.record.POS - 1
    assert region_start_index <= record_start_index

    start_offset = record_start_index - region_start_index

    region_alt = str(first_region.record.ALT[0])
    new_record_alt = region_alt[:start_offset] + str(vcf_record.ALT[0])

    vcf_record = _modify_vcf_record(vcf_record,
                                    POS=first_region.record.POS,
                                    REF=first_region.record.REF,
                                    ALT=[new_record_alt])
    return vcf_record


def _pad_vcf_record_end(vcf_record, overlap_regions: _Regions):
    """Pad end of VCF record.

    Pad ALT end with region end ALT (if any from current end backwards).
    """
    record_start_index = vcf_record.POS - 1
    record_end_index = record_start_index + len(vcf_record.REF) - 1

    last_region = overlap_regions[-1]
    overlap_end_index = record_end_index - last_region.start

    region_alt = str(last_region.record.ALT[0])
    record_alt = str(vcf_record.ALT[0])
    alt_extension = region_alt[overlap_end_index + 1:]
    new_alt = record_alt + alt_extension

    vcf_record = _modify_vcf_record(vcf_record, ALT=[new_alt])
    return vcf_record


def _find_overlap_regions(vcf_record, first_region_index, secondary_regions):
    record_start_index = vcf_record.POS - 1

    if len(vcf_record.ALT[0]) > len(vcf_record.REF):
        record_end_index = record_start_index + len(vcf_record.ALT[0]) - 1
    else:
        record_end_index = record_start_index + len(vcf_record.REF) - 1
    record_max_length = record_end_index - record_start_index + 1

    record_length_consumed = 0
    regions = iter(secondary_regions[first_region_index:])
    overlap_regions = []
    first_region_processed = False
    while record_length_consumed < record_max_length:
        region = next(regions, None)
        if region is None:
            break

        if not first_region_processed:
            start_offset = record_start_index - region.start
            record_length_consumed = len(region) - start_offset
            first_region_processed = True
        else:
            record_length_consumed += len(region)

        overlap_regions.append(region)
    return overlap_regions


def _find_start_region_index(vcf_record, secondary_regions):
    record_start_index = vcf_record.POS - 1
    region_index = bisect.bisect_left(secondary_regions, record_start_index)

    if region_index > len(secondary_regions) - 1:
        region_index = len(secondary_regions) - 1

    region = secondary_regions[region_index]
    if record_start_index < region.start:
        region_index -= 1
    return region_index


def _secondary_regions_for_base_sites(base_records) -> _Regions:
    """Identifies base sites (via index ranges) in secondary reference.
    Any secondary reference index which falls within one of these ranges overlaps an entry in the base VCF.
    """
    # Base site ranges (inclusive) in secondary indexes
    secondary_regions = []
    index_offset = 0

    for record in base_records:
        start = index_offset + record.POS - 1
        end = start + len(record.ALT[0]) - 1

        region = _Region(start, end, record=record, offset=None)
        secondary_regions.append(region)

        index_offset += len(record.ALT[0]) - len(record.REF)
    return secondary_regions


class IterPairs:
    def __init__(self, elements):
        self._elements = iter(elements)

        self._first_next_done = False
        self._first = None
        self._second = None

    def __iter__(self):
        return self

    def __next__(self):
        if self._first_next_done is False:
            self._first = None
            self._second = next(self._elements, None)

            self._first_next_done = True
            return self._safe_yield(self._first, self._second)

        self._first = self._second
        self._second = next(self._elements, None)
        return self._safe_yield(self._first, self._second)

    @staticmethod
    def _safe_yield(first, second):
        if first is None and second is None:
            raise StopIteration
        return first, second


def _secondary_regions_for_base_nonsites(secondary_site_regions: _Regions,
                                         base_records,
                                         secondary_reference_length: int) -> _Regions:
    # add offset to secondary index to get base offset
    offset = 0
    secondary_regions = []

    for record_pair, region_pair in zip(IterPairs(base_records), IterPairs(secondary_site_regions)):
        base_record = record_pair[0]
        if base_record is not None:
            offset += len(base_record.REF) - len(base_record.ALT[0])

        first_region, second_region = region_pair
        if second_region is not None and second_region.start == 0:
            continue

        at_first_region = first_region is None
        at_mid_region = first_region is not None and second_region is not None
        at_last_region = second_region is None
        start, end = None, None

        if at_first_region:
            start = 0
            end = second_region.start - 1
            offset = 0

        if at_mid_region:
            adjacent_site_regions = first_region.end + 1 == second_region.start
            if adjacent_site_regions:
                continue
            start = first_region.end + 1
            end = second_region.start - 1

        if at_last_region:
            ends_in_site_region = first_region.end + 1 == secondary_reference_length
            if ends_in_site_region:
                continue
            start = first_region.end + 1
            end = secondary_reference_length - 1

        region = _Region(start, end, record=None, offset=offset)
        secondary_regions.append(region)
    return secondary_regions


def _merge_regions(site_regions: _Regions, nonsite_regions: _Regions) -> _Regions:
    merged = []

    site_regions = iter(site_regions)
    nonsite_regions = iter(nonsite_regions)

    site_region = next(site_regions, None)
    nonsite_region = next(nonsite_regions, None)

    while site_region is not None or nonsite_region is not None:
        if site_region is not None and nonsite_region is not None:
            get_next_site = site_region.start < nonsite_region.start
        else:
            get_next_site = nonsite_region is None

        if get_next_site:
            merged.append(site_region)
            site_region = next(site_regions, None)
        else:
            merged.append(nonsite_region)
            nonsite_region = next(nonsite_regions, None)

    return merged


def _load_records(vcf_file_path):
    file_handle = open(vcf_file_path, 'r')
    vcf_reader = vcf.Reader(file_handle)
    records = list(vcf_reader)
    return records


def _dump_fasta(inferred_reference, description, file_path):
    inferred_reference = ''.join(inferred_reference)

    seq = Seq.Seq(inferred_reference, Seq.IUPAC.unambiguous_dna)
    seq_record = SeqIO.SeqRecord(seq, '', description=description)
    record = [seq_record]

    log.debug("Writing fasta reference:\n%s\n%s", file_path, description)
    SeqIO.write(record, file_path, "fasta")


def _get_inferred_reference(reference, vcf_reader):
    inferred_reference = []
    record = next(vcf_reader, None)

    start_ref_idx = -1
    end_ref_idx = -1

    for idx, x in enumerate(reference):
        within_replaced_ref = start_ref_idx < idx <= end_ref_idx
        if within_replaced_ref:
            continue

        vcf_finished = record is None
        if vcf_finished:
            inferred_reference.append(x)
            continue

        start_ref_idx = record.POS - 1
        end_ref_idx = start_ref_idx + len(record.REF) - 1

        at_allele_insert_idx = idx == start_ref_idx
        if at_allele_insert_idx:
            inferred_reference += list(str(record.ALT[0]))
            record = next(vcf_reader, None)
            continue

        inferred_reference.append(x)
    return inferred_reference


def _get_reference(prg_seq) -> typing.List:
    """Get reference by parsing a PRG sequence."""
    log.debug('Parsing PRG for reference')
    regions = prg.parse(prg_seq)

    reference = []
    for region in regions:
        region = region.alleles[0]
        reference += list(region)
    return reference
