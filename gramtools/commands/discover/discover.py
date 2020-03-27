"""
Variant discovery against inferred personalised reference.
Works on a haploid personalised reference produced by gramtools `genotype`
Any newly discovered variants get rebased, ie expressed in the original prg
    reference coordinate space.
"""
from typing import List, Dict, Iterable
import bisect
import logging
import collections
import json

from pysam import VariantFile, VariantRecord
import cortex
from Bio import SeqIO

from gramtools.commands.paths import DiscoverPaths
from .region_mapper import RegionMapper, ChromSizes, _Regions

log = logging.getLogger("gramtools")


def run(args):
    log.info("Start process: discover")
    disco_paths = DiscoverPaths(args.disco_dir, args.geno_dir, args.force)
    disco_paths.setup()

    enforce_genotyping_was_haploid(disco_paths)

    # call variants using `cortex`
    log.debug("Running cortex")
    cortex.calls(
        disco_paths.pers_ref, disco_paths.reads_files, disco_paths.discov_vcf_cortex
    )

    chrom_sizes: ChromSizes = dict()
    with open(disco_paths.pers_ref) as genome_file:
        for seq_record in SeqIO.parse(genome_file, "fasta"):
            chrom_sizes[seq_record.id] = len(seq_record.seq)

    #  Convert coordinates in personalised reference space to coordinates in (original) prg space.
    log.debug("Rebasing vcf")
    rebased_vcf_records = _rebase_vcf(disco_paths, chrom_sizes)

    if rebased_vcf_records is None:
        log.debug("Rebased VCF does not contain records")
        log.info("End process: discover")
        return

    _dump_rebased_vcf(rebased_vcf_records, disco_paths)

    log.info(
        "End process: discover. " "Rebased vcf in {}".format(disco_paths.final_vcf)
    )


def _rebase_vcf(disco_paths: DiscoverPaths, chrom_sizes, check_records=True):
    """Rebase a vcf so that it uses same reference as base_vcf.
    Input:
     discovery.vcf                   personalised_ref.vcf
      |                               |
     personalised_ref.fasta          base_ref.fasta

    Output:
     discovery.vcf
      |
     base_ref.fasta
    """
    base_vcf_file_path = disco_paths.geno_vcf

    if check_records:
        var_unplaced_records = []
        inferred_refs = load_multi_fasta(disco_paths.pers_ref)

    base_records = VariantFile(base_vcf_file_path).fetch()
    derived_records = VariantFile(disco_paths.discov_vcf_cortex).fetch()

    region_mapper = RegionMapper(base_records, chrom_sizes)
    flagged_personalised_regions = region_mapper.get_mapped()

    new_vcf_records = []
    for vcf_record in derived_records:
        chrom_key = vcf_record.chrom

        if check_records:
            if (
                check_ref_consistent(
                    vcf_record, inferred_refs[chrom_key], var_unplaced_records
                )
                is False
            ):
                continue  # Do not process inconsistent records

        regions_list = flagged_personalised_regions.get(chrom_key, None)
        assert regions_list is not None, (
            f"Error: ref ID {chrom_key} in vcf record {vcf_record} is not"
            f" present in vcf from infer {base_vcf_file_path}"
        )

        new_record = _rebase_vcf_record(vcf_record, regions_list)

    if check_records and len(var_unplaced_records) > 0:
        log.warning(
            "The following new variant records were skipped: {} \n"
            "Reasons: pos is not a valid position in the inferred reference sequence,"
            "or ref does not coincide with inferred reference sequence.".format(
                "\t".join(var_unplaced_records)
            )
        )

    return new_vcf_records


def _dump_rebased_vcf(records: List[VariantRecord], disco_paths: DiscoverPaths):
    template_vcf = VariantFile(disco_paths.discov_vcf_cortex)
    if list(template_vcf.header.contigs) == list():  # Need to add contig headers
        genotyped_vcf = VariantFile(disco_paths.geno_vcf)
        for contig in genotyped_vcf.header.contigs:
            template_vcf.header.add_line(f"##contig=<ID={contig}>")

    output_vcf = VariantFile(disco_paths.final_vcf, "w", header=template_vcf.header)
    for record in records:
        output_vcf.write(record)


def _modify_vcf_record(vcf_record, **new_attributes):
    for attribute, value in new_attributes.items():
        if getattr(vcf_record, attribute, None):
            setattr(vcf_record, attribute, value)

    return vcf_record


def _rebase_vcf_record(
    vcf_record: VariantRecord, personalised_reference_regions: _Regions
):
    """Change `vcf_record` to be expressed relative to a different reference."""

    # Get index of region containing the start of the `vcf_record`.
    region_index = _find_start_region_index(vcf_record, personalised_reference_regions)

    consumed_reference = 0  # Position in the inferred ref. sequence
    reference_length = len(vcf_record.ref)
    ref_seq_left = True

    # we will build the rebased_ref as we traverse the regions. We will use the vcf_record alt and pre-pend and/or post-pend to it if necessary.
    rebased_ref, rebased_alt = "", str(vcf_record.alts[0])

    # Let's rebase the position straight away
    first_region = personalised_reference_regions[region_index]

    # Case: hitting variant region. We rebase at the beginning of the variant region.
    if first_region.is_site:
        rebased_pos = first_region.base_pos

        # We also straight away pre-pend any preceding variation relative to the base ref
        if vcf_record.pos > first_region.inf_pos:
            record_inset = vcf_record.pos - first_region.inf_pos
            rebased_alt = first_region.vcf_record_alt[:record_inset] + rebased_alt

    # Case: hitting non-variant region. We rebase at where the vcf_record starts, in base ref coordinates.
    else:
        rebased_pos = first_region.base_pos + (vcf_record.pos - first_region.inf_pos)

    while ref_seq_left:
        region = personalised_reference_regions[region_index]
        # Check how much of the vcf_record ref (inferred reference) can be consumed by the current region.
        # If the current region can consume at least what is left of the vcf_record ref, loop ends.
        # NOTE that region.length is 'overloaded': if a non-var region, it is the fixed interval between var regions.
        # If a var region, it is the inferred_vcf record's alt length (ref and alt lengths can differ).
        consumable = region.length - (
            vcf_record.pos + consumed_reference - region.inf_pos
        )

        if consumable >= (reference_length - consumed_reference):
            ref_seq_left = False
            to_consume = reference_length - consumed_reference
        else:
            to_consume = consumable

        if region.is_site:
            rebased_ref += region.vcf_record_ref

        else:
            # We can use the vcf_record's ref, as that is also the base ref sequence- because we are in a non-variant site.
            rebased_ref += vcf_record.ref[
                consumed_reference : consumed_reference + to_consume
            ]

        consumed_reference += to_consume
        region_index += 1

    assert consumed_reference == len(vcf_record.ref)

    # Deal with the last region: post-pend any sequence in alt record if we finish in a variant site.
    if region.is_site:
        cur_pos = vcf_record.pos + consumed_reference
        # The inset will be < 0 if there is a part of the (inferred vcf record's) alt which has not been
        inset = cur_pos - (region.inf_pos + region.length)
        if inset < 0:
            rebased_alt += region.vcf_record_alt[inset:]

    vcf_record = _modify_vcf_record(
        vcf_record, pos=rebased_pos, ref=rebased_ref, alts=[rebased_alt]
    )

    return vcf_record


## Return the marked `_Region` containing the position referred to by `vcf_record`.
def _find_start_region_index(vcf_record: VariantRecord, marked_regions):
    record_start_pos = vcf_record.pos

    # Binary search on a sorted array. '<' comparator is overloaded in `_Region` definition to allow the search.
    # The produced index is:
    # * If a region with start pos == record_start_pos is found, returns that index (specifically, the first occurrence)
    # * If such a region is not found, returns the index of the region with start pos **strictly** greater than record_start_pos
    # IN THAT CASE, we need to return the region index one smaller; so that the region's interval includes record_start_pos.
    region_index = bisect.bisect_left(marked_regions, record_start_pos)

    #  Case: `record_start_index` is larger than any _Region start.
    if region_index > len(marked_regions) - 1:
        region_index = len(marked_regions) - 1

    selected_region = marked_regions[region_index]

    # Case: record_start_pos was not found exactly.
    if selected_region.inf_pos > record_start_pos:
        region_index -= 1

    return region_index


def enforce_genotyping_was_haploid(disco_paths: DiscoverPaths):
    report_path = disco_paths.geno_report
    with open(report_path) as f_in:
        genotype_report = json.load(f_in)
    gt_ploidy = genotype_report["ploidy"]
    if gt_ploidy != "haploid":
        log.error(
            f"Discover currently supports haploid genotyping only.\n"
            f"The report from genotype at {report_path} states genotyping ran in {gt_ploidy} mode.\n"
            f"If your organism really is {gt_ploidy}, discover variants from haploid,"
            f"then run genotype in {gt_ploidy} mode :-)"
        )
        exit(1)


def check_ref_consistent(
    vcf_record: VariantRecord, inferred_sequence: str, var_unplaced_record: List[str]
):
    """Check the variant call is properly placed by the variant caller"""
    position = vcf_record.pos
    if (
        len(inferred_sequence) < position
        or vcf_record.ref
        != inferred_sequence[position - 1 : position - 1 + len(vcf_record.ref)]
    ):
        var_unplaced_records.append(str(vcf_record))
        return False
    return True


## Generator to always yield the first allele (ref) of a variant site.
def _alwaysRef():
    while True:
        yield 0


def load_multi_fasta(fasta_fname):
    # Load the record fully in memory
    with open(fasta_fname) as r:
        records = collections.OrderedDict()
        record = ""
        new_recomb_name = ""

        for line in r:
            if line[0] == ">":
                if new_recomb_name != "":
                    records[new_recomb_name] = record
                new_recomb_name = line[1:].strip().split()[0]
                record = ""

            else:
                record += line.replace("\n", "")

        records[new_recomb_name] = record

    return records
