## @file
# Variant discovery against inferred personalised reference.
# The use case is the following: a set of reads is mapped to a prg using gram::commands::quasimap.
# Then a haploid personalised reference is inferred using `infer.py`.
# Finally the same set of reads can be mapped to this personalised reference for variant discovery using this module.
# The output vcf file gets rebased to express variants in the original prg (rather than personalised reference) coordinate space.

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
from .. import prg_regions_parser

log = logging.getLogger('gramtools')


def parse_args(common_parser, subparsers):
    parser = subparsers.add_parser('discover',
                                   parents=[common_parser])
    parser.add_argument('--gram-directory',
                        help='',
                        type=str,
                        required=True)
    parser.add_argument('--inferred-vcf',
                        help='The vcf file corresponding to the inferred path through the PRG.',
                        type=str,
                        required=True)
    parser.add_argument('--reads',
                        help='Reads file for variant discovery with respect to the inferred personalised reference.',
                        type=str,
                        required=True)
    parser.add_argument('--output-vcf',
                        help='File path for variant calls.',
                        type=str,
                        required=False)


def run(args):
    log.info('Start process: discover')
    _paths = paths.generate_discover_paths(args)
    os.mkdir(_paths['tmp_directory'])

    with open(_paths['prg'], 'r') as file_handle:
        prg_seq = file_handle.read()
    # `reference` is all non-variant parts of prg + first allele of each variant site of prg
    reference = _get_reference(prg_seq)
    fasta_description = "PRG base reference for gramtools"
    _dump_fasta(reference, fasta_description, _paths['base_reference'])

    file_handle = open(args.inferred_vcf, 'r')
    vcf_reader = vcf.Reader(file_handle)

    inferred_reference = _get_inferred_reference(reference, vcf_reader)
    fasta_description = "inferred personal reference generated using gramtools"
    _dump_fasta(inferred_reference, fasta_description, _paths['inferred_reference'])

    # call variants using `cortex`
    cortex.calls(_paths['inferred_reference'], args.reads, _paths['cortex_vcf'])

    inferred_reference_length = len(inferred_reference)
    # Convert coordinates in personalised reference space to coordinates in (original) prg space.
    rebased_vcf_records = _rebase_vcf(args.inferred_vcf,
                                      inferred_reference_length,
                                      _paths['cortex_vcf'])
    if rebased_vcf_records is None:
        log.debug("Rebased VCF does not contain records")
        shutil.rmtree(_paths['tmp_directory'])
        log.info('End process: discover')
        return

    template_vcf_file_path = _paths['cortex_vcf']
    _dump_rebased_vcf(rebased_vcf_records, _paths['rebase_vcf'], template_vcf_file_path)

    vcf_files = [_paths['rebase_vcf']]
    # Cluster together variants in close proximity.
    cluster = cluster_vcf_records.vcf_clusterer.VcfClusterer(vcf_files, _paths['base_reference'], args.output_vcf)
    cluster.run()

    shutil.rmtree(_paths['tmp_directory'])
    log.info('End process: discover')

## Rebase a vcf so that it uses same reference as base_vcf.
#
# Here is the use case. We have:
#  new_vcf                         new_reference_vcf
#   |                               |
#  new_reference                   old_reference
#
# And we want:
#  new_vcf
#   |
#  old_reference
#
# In the context of gramtools:
# * `new_vcf` is the vcf produced by cortex, ie variant calls against personalised reference.
# * `new_reference_vcf` is the vcf describing the personalised reference with respect to the original prg.
# * `old_reference` is the original prg. The first allele of each variant site in the prg is taken for making this.
# * `new_reference` is the personalised reference; it is the original prg with inferred variant alleles instead.
#
# Thus what rebasing achieves is to express the variant calls in terms of the original prg coordinates.
def _rebase_vcf(base_vcf_file_path,
                personalised_reference_length: int,
                personalised_vcf_file_path):

    try:
        # Load the variant calls with respect to the personalised reference.
        variant_call_records = _load_records(personalised_vcf_file_path)
    except StopIteration:
        log.warning("Cortex VCF does not contain novel variants")
        return None
    base_records = _load_records(base_vcf_file_path)

    # Given a variant call position w.r.t the personalised reference, we will need to know how this position maps to the original prg.
    flagged_personalised_regions = _flag_personalised_reference_regions(base_records, personalised_reference_length)

    new_vcf_records = [_rebase_vcf_record(vcf_record, flagged_personalised_regions)
                       for vcf_record in variant_call_records]
    return new_vcf_records


def _dump_rebased_vcf(records, rebase_file_path, template_vcf_file_path):
    template_handle = open(template_vcf_file_path, 'r')
    template_vcf = vcf.Reader(template_handle)

    rebase_handle = open(rebase_file_path, 'w')
    writer = vcf.Writer(rebase_handle, template_vcf)
    for record in records:
        writer.write_record(record)
    writer.close()

## Class that marks a region in a reference (call it 2) with a relationship to another reference (call it 1).
# By 'earlier' we imply that reference 2 is based off of reference 1.
# A `_Region` either marks a region in ref 2 within a variant region in ref 1, or a region in ref 2 within an invariant region in ref 1.
class _Region:
    def __init__(self, start, end, offset, record):
        ## Start coordinate in personalised reference.
        self.start = start
        ## End coordinate in personalised reference.
        self.end = end
        ## Holds the total difference in coordinates between ref 2 and ref 1. Only populated if region marks a non-variant site in ref 1.
        self.offset = offset
        ## Holds vcf record for this region of ref 2 relative to ref 1. Only populated if region marks a variant site in ref 1.
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



## Marks the personalised reference with `Regions` flagging a relationship with the base reference.
# There are two types of `Regions`: those inside variant sites in the prg, and those outside.
def _flag_personalised_reference_regions(base_records, secondary_reference_length) -> _Regions:
    site_regions = mark_base_site_regions(base_records)
    nonsite_regions = mark_base_nonsite_regions(site_regions,
                                                base_records,
                                                secondary_reference_length)

    all_personalised_ref_regions = _merge_regions(site_regions, nonsite_regions)
    return all_personalised_ref_regions


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

## Change `vcf_record` to be expressed relative to a different reference.
# @param vcf_record a representation of a vcf entry;
# @param personalised_reference_regions a set of `_Regions` marking the personalised reference and its relationship to the prg reference.
def _rebase_vcf_record(vcf_record, personalised_reference_regions: _Regions):
    # Get index of region containing the start of the `vcf_record`.
    first_region_index = _find_start_region_index(vcf_record, personalised_reference_regions)
    overlap_regions = _find_overlap_regions(vcf_record, first_region_index, personalised_reference_regions)

    # TODO: what about start and end in nonsites and sites in the middle? below condition probably not general enough
    overlap_sites = any(region.is_site for region in overlap_regions)

    # Case: the variant region does not overlap any variant sites in the prg reference.
    # All we need to do is modify the coordinates of the record to map to the prg reference.
    if not overlap_sites:
        offset_sum = sum(region.offset for region in overlap_regions) # TODO: if we do not overlap any sites, we should overlap only ONE site?
        new_pos = vcf_record.POS + offset_sum

        vcf_record = _modify_vcf_record(vcf_record, POS=new_pos)
        return vcf_record

    if overlap_regions[0].is_site:
        vcf_record = _pad_vcf_record_start(vcf_record, overlap_regions)

    if overlap_regions[-1].is_site:
        vcf_record = _pad_vcf_record_end(vcf_record, overlap_regions)

    return vcf_record


## Pad start of VCF record with what the original variant has.
# The POS and REf of `vcf_record` get set to what they are in the original prg.
def _pad_vcf_record_start(vcf_record, overlap_regions: _Regions):
    original_variant_region = overlap_regions[0]

    record_start_index = vcf_record.POS - 1
    original_region_start_index = original_variant_region.record.POS - 1
    assert original_region_start_index <= record_start_index

    start_offset = record_start_index - original_region_start_index

    region_alt = str(original_variant_region.record.ALT[0])
    new_record_alt = region_alt[:start_offset] + str(vcf_record.ALT[0])

    # REF gets set to what the original (base; prg-level) REF is, as stored in the `_Region` object.
    vcf_record = _modify_vcf_record(vcf_record,
                                    POS=original_variant_region.record.POS,
                                    REF=original_variant_region.record.REF,
                                    ALT=[new_record_alt])
    return vcf_record


## Pad end of VCF record with what the original variant has.
def _pad_vcf_record_end(vcf_record, overlap_regions: _Regions):
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

## Find all `_Region`s overlapped by a `vcf_record`.
# @return the overlapped `_Region`s.
def _find_overlap_regions(vcf_record, first_region_index, marked_regions):
    record_start_index = vcf_record.POS - 1

    # Take the largest region alluded to in the `vcf_record`. This region is used to query which `_Region`s are overlapped.
    # TODO: should it not be the REF length that is used as record length, rather than max of REF and ALT?
    if len(vcf_record.ALT[0]) > len(vcf_record.REF):
        record_end_index = record_start_index + len(vcf_record.ALT[0]) - 1
    else:
        record_end_index = record_start_index + len(vcf_record.REF) - 1
    record_max_length = record_end_index - record_start_index + 1

    record_length_consumed = 0
    regions = iter(marked_regions[first_region_index:])
    overlap_regions = []
    process_first_region = True

    # Keep appending regions while we have not traversed the personalised reference with the full length of the record.
    while record_length_consumed < record_max_length:
        region = next(regions, None)
        if region is None:
            break

        # When processing the first region, the `record_start_index` can be greater than the start of the region.
        # We do not consume the whole region thus, but only the region between the `record_start_index` and the region's end.
        if process_first_region:
            start_offset = record_start_index - region.start
            record_length_consumed = len(region) - start_offset
            process_first_region = False
        else:
            record_length_consumed += len(region)

        overlap_regions.append(region)
    return overlap_regions

## Return the marked `_Region` containing the position referred to by `vcf_record`.
def _find_start_region_index(vcf_record, marked_regions):
    record_start_index = vcf_record.POS - 1
    # Bisect-left is a binary search on a sorted array. The '<' comparator is overloaded in `_Region` definition to allow the search.
    region_index = bisect.bisect_left(marked_regions, record_start_index)

    # Case: We have picked out the very last region.
    if region_index > len(marked_regions) - 1:
        region_index = len(marked_regions) - 1

    region = marked_regions[region_index]
    # This part is important. The `bisect_left` routine can produce the index of the region **just to the right**
    # of the region containing the `vcf_record`. This will in fact always be the case when the vcf_record position is
    # after the region's start. Then we decrement the index.
    if record_start_index < region.start:
        region_index -= 1
    return region_index


## Flags `_Regions` in the personalised reference which are in variant sites in the original PRG.
# Any personalised reference position which falls within one of these `Regions` overlaps an entry in the original VCF.
# Thus, the `_Regions` are marked with the base vcf record.
def mark_base_site_regions(base_records) -> _Regions:
    # Base site ranges (inclusive) in personalised reference
    secondary_regions = []
    index_offset = 0

    for record in base_records:
        start = index_offset + record.POS - 1
        end = start + len(record.ALT[0]) - 1

        region = _Region(start, end, record=record, offset=None)
        secondary_regions.append(region)

        # The offset keeps track of how much larger or smaller the personalised reference is relative to the prg reference.
        index_offset += len(record.ALT[0]) - len(record.REF)
    return secondary_regions


## A class for iterating through elements in pairs.
# There are two edge cases:
# * The very first iteration sets up a pair of an empty element and the first element of `elements`.
# * The very last iteration sets up a pair of the very last element of `elements` and an empty element.
# @param elements the elements to iterate through in pairs.
class IterPairs:
    def __init__(self, elements):
        self._elements = iter(elements) # Sets up iterator on `elements` so we can call next() on them.
        ## This variable flags whether we have iterated past the very first pair.
        self._first_next_done = False
        self._first = None
        self._second = None

    def __iter__(self):
        return self

    def __next__(self):
        if self._first_next_done is False: # Iterate through very first pair: returns the very first element paired to an empty element.
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
            raise StopIteration # Nothing left to give, stop there
        return first, second


## Flags `_Regions` in the personalised reference which are in non-variant sites in the original PRG.
# Use `IterPairs` class because we need the information on two consecutive site region to compute `_Region` coordinates.
def mark_base_nonsite_regions(site_regions: _Regions,
                              base_records,
                              secondary_reference_length: int) -> _Regions:
    # add offset to secondary index to get base offset
    offset = 0
    secondary_regions = []

    for record_pair, region_pair in zip(IterPairs(base_records), IterPairs(site_regions)):
        base_record = record_pair[0]
        if base_record is not None: # This is only false once, when we encounter the very first record. Then offset must stay 0.
            offset += len(base_record.REF) - len(base_record.ALT[0])

        first_region, second_region = region_pair # Unpack the two regions

        # If the first position of the personalised reference corresponds to a variant site in the prg, there is nothing to record.
        if second_region is not None and second_region.start == 0:
            continue

        at_first_region = first_region is None # `IterPairs` produces an empty `first_region` for the very first pair.
        at_mid_region = first_region is not None and second_region is not None
        at_last_region = second_region is None # `IterPairs` produces an empty `second_region` for the very last pair.
        start, end = None, None

        if at_first_region:
            start = 0
            end = second_region.start - 1
            offset = 0

        if at_mid_region:
            adjacent_site_regions = first_region.end + 1 == second_region.start
            # If there is no gap between the pair of site regions, there is nothing to record.
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

## Merges the `Regions` (site and nonsite) by ordering them using position coordinates.
def _merge_regions(site_regions: _Regions, nonsite_regions: _Regions) -> _Regions:
    merged = []

    site_regions = iter(site_regions)
    nonsite_regions = iter(nonsite_regions)

    site_region = next(site_regions, None)
    nonsite_region = next(nonsite_regions, None)

    while site_region is not None or nonsite_region is not None:
        if site_region is not None and nonsite_region is not None:
            # Case 1: we have both a site_region and a nonsite_region.
            # `get_next_site` flags whether we next need to deal with the site region (`true`) or the nonsite region (`false`).
            get_next_site = site_region.start < nonsite_region.start
        else:
            # Case 2: we have only one type of `Region` left.
            # `get_next_site` flags whether that `Region` is a site region (`true`) or the nonsite region.
            get_next_site = nonsite_region is None

        if get_next_site: # Append the site_region, and get the next one.
            merged.append(site_region)
            site_region = next(site_regions, None)
        else:
            merged.append(nonsite_region)
            nonsite_region = next(nonsite_regions, None)

    return merged

## Produce a list of vcf records parsed from `vcf_file_path`.
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

## Compute a new reference from an old reference and a set of vcf records.
# @param reference the old reference.
# @param vcf_reader set of vcf records describing variants against the old reference.
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
    regions = prg_regions_parser.parse(prg_seq)

    reference = []
    for region in regions:
        region = region.alleles[0]
        reference += list(region)
    return reference
