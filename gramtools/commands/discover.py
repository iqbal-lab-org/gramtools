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
from .. import prg_local_parser

REF_READ_CHUNK_SIZE = 1000000
DNA = {'A', 'C', 'G', 'T'}

log = logging.getLogger('gramtools')


def parse_args(common_parser, subparsers):
    parser = subparsers.add_parser('discover',
                                   parents=[common_parser])
    parser.add_argument('--gram-dir','--gram-directory',
                        help='',
                        type=str,
                        dest='gram_directory',
                        required=True)
    parser.add_argument('--infer-dir',
                        help='The directory containing the outputs from `infer` command',
                        type=str,
                        required=True)
    parser.add_argument('--reads',
                        help='Reads file for variant discovery with respect to the inferred personalised reference.',
                        type=str,
                        required=True)
    parser.add_argument('--discover-dir',
                        help='Directory to output results from `discover` command'
                             'Defaults to inside "gram-dir"',
                        type=str,
                        required=False)


def run(args):
    log.info('Start process: discover')
    _paths = paths.generate_discover_paths(args)
    # os.mkdir(_paths['tmp_directory'])


    # call variants using `cortex`
    cortex.calls(_paths['inferred_fasta'], args.reads, _paths['cortex_vcf'])


    # Get the length of the inferred fasta reference; was computed at infer stage.
    with open(_paths["inferred_ref_size"]) as f:
        inferred_reference_length = int(f.read())

    # Convert coordinates in personalised reference space to coordinates in (original) prg space.
    rebased_vcf_records = _rebase_vcf(_paths["inferred_vcf"],
                                      inferred_reference_length,
                                      _paths['cortex_vcf'])
    if rebased_vcf_records is None:
        log.debug("Rebased VCF does not contain records")
        shutil.rmtree(_paths['tmp_directory'])
        log.info('End process: discover')
        return

    template_vcf_file_path = _paths['cortex_vcf']
    _dump_rebased_vcf(rebased_vcf_records, _paths['rebase_vcf'], template_vcf_file_path)

    vcf_files = [_paths['rebased_vcf']]
    # Cluster together variants in close proximity.
    cluster = cluster_vcf_records.vcf_clusterer.VcfClusterer(vcf_files, _paths['base_reference'], _paths['final_output_vcf'])
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

## Class that allows mapping between a vcf record in one REF coordinate to a vcf record in another REF coordinate system.
# Call the two references coordinates ref 2 and ref 1.
# A `_Region` either marks a region in ref 2 within a variant region in ref 1, or a region in ref 2 within an invariant region in ref 1.
class _Region:
    def __init__(self, base_POS, inf_POS, length, vcf_record_REF = None, vcf_record_ALT = None):
        ## Start coordinate in base reference.
        self.base_POS = base_POS
        ## Start coordinate in inferred reference.
        self.inf_POS = inf_POS
        ## Holds the length of the region
        self.length = length

        ## If in a variant site, the only parts of the vcf_record that we need to track are the REF and ALT sequences.
        self.vcf_record_REF = vcf_record_REF

        self.vcf_record_ALT = vcf_record_ALT

    @property
    def is_site(self):
        return self.vcf_record_REF is not None

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __repr__(self):
        return str(self.__dict__)

    def __lt__(self, other):
        if isinstance(other, _Region):
            return self.inf_POS < other.inf_POS
        else:
            return self.inf_POS < other

    def __len__(self):
        return self.end - self.start + 1


_Regions = typing.List[_Region]



## Marks the personalised reference with `Regions` flagging a relationship with the base reference.
# There are two types of `Regions`: those inside variant sites in the prg, and those outside.
def _flag_personalised_reference_regions(base_records, secondary_reference_length) -> _Regions:

    all_personalised_ref_regions = []

    base_ref_pos = 1
    inf_ref_pos = 1

    for vcf_record in base_records:
        # Build non-variant region
        if vcf_record.POS > base_ref_pos:
            non_var_region = _Region(base_POS = base_ref_pos, inf_POS = inf_ref_pos,
                                     length = vcf_record.POS - base_ref_pos)
            all_personalised_ref_regions.append(non_var_region)
            base_ref_pos += non_var_region.length
            inf_ref_pos += non_var_region.length

        # Build variant region
        var_region = _Region(base_POS = base_ref_pos, inf_POS = inf_ref_pos,
                             length = len(vcf_record.ALT[0]), vcf_record_REF = vcf_record.REF,
                             vcf_record_ALT = str(vcf_record.ALT[0]))
        all_personalised_ref_regions.append(var_region)

        base_ref_pos += len(vcf_record.REF)
        inf_ref_pos += var_region.length

    # End game: deal with the last non-var region if there is one.
    # We test `inf_ref_pos` because we have the inferred reference length (= secondary reference)
    if inf_ref_pos <= secondary_reference_length:
        non_var_region = _Region(base_POS = base_ref_pos, inf_POS = inf_ref_pos,
                                 length = secondary_reference_length - inf_ref_pos + 1)
        all_personalised_ref_regions.append(non_var_region)

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
# The key idea is to progress through the REF sequence of the vcf_record to rebase, mapping it to the variant/non-variant space of the personalised reference.
def _rebase_vcf_record(vcf_record, personalised_reference_regions: _Regions):

    # Get index of region containing the start of the `vcf_record`.
    region_index = _find_start_region_index(vcf_record, personalised_reference_regions)


    consumed_reference = 0  # Position in the inferred ref. sequence
    reference_length = len(vcf_record.REF)

    # We will iterate over regions until we have consumed all of the vcf_record REF sequence.
    journey_continues = True

    # we will build the rebased_REF as we traverse the regions. We will use the vcf_record ALT and pre-pend and/or post-pend to it if necessary.
    rebased_REF, rebased_ALT = "", str(vcf_record.ALT[0])


    # Let's rebase the position straight away
    first_region = personalised_reference_regions[region_index]

    # Case: hitting variant region. We rebase at the beginning of the variant region.
    if first_region.is_site:
        rebased_POS = first_region.base_POS

        # We also straight away pre-pend any preceding variation relative to the base REF
        if vcf_record.POS > first_region.inf_POS:
            record_inset = vcf_record.POS - first_region.inf_POS
            rebased_ALT = first_region.vcf_record_ALT[:record_inset] + rebased_ALT

    # Case: hitting non-variant region. We rebase at where the vcf_record starts, in base ref coordinates.
    else:
        rebased_POS = first_region.base_POS + (vcf_record.POS - first_region.inf_POS)


    while journey_continues:
        region = personalised_reference_regions[region_index]

        # This is the key step. We check how much of the vcf_record REF (inferred reference) can be consumed by the current region.
        # If the current region can consume at least what is left of the vcf_record REF, loop ends.
        # NOTE that region.length is 'overloaded': if a non-var region, it is the fixed interval between var regions.
        # If a var region, it is the inferred_vcf record's ALT length (REF and ALT lengths can differ).
        consumable = region.length - (vcf_record.POS + consumed_reference - region.inf_POS)  # The amount of the inferred reference that can be consumed by the region

        if consumable >= (reference_length - consumed_reference):
            journey_continues = False
            to_consume = reference_length - consumed_reference
        else:
            to_consume = consumable


        if region.is_site:
            rebased_REF += region.vcf_record_REF   # Copy the whole base REF for the variant site

        else:
            # We can use the vcf_record's REF, as that is also the base REF sequence- because we are in a non-variant site.
            rebased_REF += vcf_record.REF[consumed_reference:consumed_reference + to_consume]

        consumed_reference += to_consume  # Log the progression

        region_index += 1

    # Sanity check: we have consumed the whole of vcf_record.REF.
    assert consumed_reference == len(vcf_record.REF)

    # End game: Deal with the last region. We need to post-pend any succeeding sequence in the ALT record if we finish in a variant site.
    if region.is_site:
        cur_pos = vcf_record.POS + consumed_reference
        # The inset will be < 0 if there is a part of the (inferred vcf record's) ALT which has not been
        inset = cur_pos - (region.inf_POS + region.length)
        if inset < 0:
            rebased_ALT += region.vcf_record_ALT[inset:]

    vcf_record = _modify_vcf_record(vcf_record, POS = rebased_POS, REF = rebased_REF, ALT = [rebased_ALT])

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
    record_start_POS = vcf_record.POS

    # Binary search on a sorted array. The '<' comparator is overloaded in `_Region` definition to allow the search.
    # The produced index is:
    # * If a region with start POS == record_start_POS is found, returns that index (specifically, the first occurrence)
    # * If such a region is not found, returns the index of the region with start POS **strictly** greater than record_start_POS
    # IN THAT CASE, we need to return the region index one smaller; so that the region's interval includes record_start_POS.
    region_index = bisect.bisect_left(marked_regions, record_start_POS)

    # Case: `record_start_index` is larger than any _Region start.
    if region_index > len(marked_regions) - 1:
        region_index = len(marked_regions) - 1

    selected_region = marked_regions[region_index]

    # Case: record_start_POS was not found exactly.
    if selected_region.inf_POS > record_start_POS:
        region_index -= 1

    return region_index



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


## Compute a new reference from an base reference and a set of vcf records.
# @param reference the 'base' reference: non-variant sites of the prg + first allele of each variant site.
# @param vcf_reader set of vcf records describing variants inferred against the old reference.
def _produce_inferred_reference(reference_fpath, vcf_reader, fasta_fpath, description):
    inferred_reference_length = 0

    reference = _refParser(reference_fpath)
    inferred_reference = prg_local_parser.FastaWriter(fasta_fpath, description)

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
            inferred_reference_length += 1
            continue

        # This update only actually occurs once we have passed a record
        start_ref_idx = record.POS - 1
        end_ref_idx = start_ref_idx + len(record.REF) - 1

        at_allele_insert_idx = idx == start_ref_idx

        # Let's insert the genotyped allele
        if at_allele_insert_idx:
            sample = record.samples[0]
            genotype = sample.gt_alleles

            # Case: no genotype. Let's pick REF allele.
            if set(genotype) == {None}:
                allele_index = 0

            # Case: single haploid genotype. Pick that allele.
            # OR Case: there are >1 calls. Just pick the first.
            else:
                allele_index = int(genotype[0])


            if allele_index == 0:
                template = str(record.REF)
            else:
                template = str(record.ALT[allele_index - 1])

            for alt_char in template:
                inferred_reference.append(alt_char)

            inferred_reference_length += len(template)

            record = next(vcf_reader, None)
            continue

        # Add the non-variant regions
        inferred_reference.append(x)
        inferred_reference_length += 1

    inferred_reference.close()

    return inferred_reference_length


## Generator to always yield the first allele (REF) of a variant site.
def _alwaysRef():
    while True:
        yield 0


## Generator to parse fasta ref sequence in blocks.
def _refParser(ref_fpath, chunk_size = REF_READ_CHUNK_SIZE):
    with open(ref_fpath) as fhandle:
        next(fhandle)  # Skip header line
        while True:
            chars = fhandle.read(chunk_size)
            if len(chars) == 0:
                return

            for char in chars:
                if char.upper() not in DNA:
                    continue
                yield char


