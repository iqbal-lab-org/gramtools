"""
Variant discovery against inferred personalised reference.
Works on a haploid personalised reference produced by gramtools `genotype`
Any newly discovered variants get rebased, ie expressed in the original prg
    reference coordinate space.
"""
from typing import List
import logging
import json

from pysam import VariantFile, VariantRecord
import cortex.calls as cortex

from gramtools.commands.paths import DiscoverPaths
from gramtools.commands.common import load_fasta
from .seq_region_map import (
    Chrom,
    ChromSizes,
    SeqRegionsMap,
    SeqRegionMapper,
    SearchableSeqRegionsMap,
    BisectTarget,
)

log = logging.getLogger("gramtools")


def run(args):
    log.info("Start process: discover")
    disco_paths = DiscoverPaths(args.disco_dir, args.geno_dir, args.force)
    disco_paths.setup()

    enforce_genotyping_was_haploid(disco_paths)

    # call variants using `cortex`
    log.debug("Running cortex")
    cortex.run(
        disco_paths.pers_ref, disco_paths.reads_files, disco_paths.discov_vcf_cortex
    )

    chrom_sizes: ChromSizes = load_fasta(disco_paths.pers_ref, sizes_only=True)

    #  Convert coordinates in personalised reference space to coordinates in (original) prg space.
    log.debug("Rebasing vcf")
    rebased_vcf_records = _rebase_vcf(disco_paths, chrom_sizes)
    num_records = len(rebased_vcf_records)

    _dump_rebased_vcf(rebased_vcf_records, disco_paths)

    log.info(f"Found {num_records} variants. " f"Final vcf in {disco_paths.final_vcf}")
    log.info("End process: discover.")


def _rebase_vcf(disco_paths: DiscoverPaths, chrom_sizes, check_records=True):
    """Rebase a vcf so that it uses same reference as base_vcf.
    (* for not an input/output, just for illustration)
    Input:
     discovery.vcf                   personalised_ref.vcf
      |                               |
     personalised_ref.fasta          *base_ref.fasta

    Output:
     discovery.vcf
      |
     *base_ref.fasta
    """
    if check_records:
        var_unplaced_records = []
        inferred_refs = load_fasta(disco_paths.pers_ref)

    _add_contig_lines(disco_paths)
    base_records = VariantFile(disco_paths.geno_vcf).fetch()
    derived_records = VariantFile(disco_paths.discov_vcf_cortex).fetch()

    region_map: SeqRegionsMap = SeqRegionMapper(base_records, chrom_sizes).get_map()
    region_searcher = SearchableSeqRegionsMap(region_map)

    new_vcf_records = []
    for vcf_record in derived_records:
        chrom_key = vcf_record.chrom

        if check_records:
            if not check_ref_consistent(
                vcf_record, inferred_refs[chrom_key], var_unplaced_records
            ):
                continue  # Do not process inconsistent records

        new_vcf_records.append(
            _rebase_vcf_record(vcf_record, chrom_key, region_searcher)
        )

    if check_records and len(var_unplaced_records) > 0:
        log.warning(
            f"{len(var_unplaced_records)} new variant records were skipped, "
            f"because record pos and ref do not coincide with personalised reference"
        )
        log.debug("Skipped records: {}".format("\n".join(var_unplaced_records)))

    return new_vcf_records


def _dump_rebased_vcf(records: List[VariantRecord], disco_paths: DiscoverPaths):
    template_vcf = VariantFile(disco_paths.discov_vcf_cortex)
    output_vcf = VariantFile(disco_paths.final_vcf, "w", header=template_vcf.header)
    for record in records:
        output_vcf.write(record)


def _modify_vcf_record(vcf_record, **new_attributes):
    for attribute, value in new_attributes.items():
        if getattr(vcf_record, attribute, None):
            setattr(vcf_record, attribute, value)

    return vcf_record


def _rebase_vcf_record(
    vcf_record: VariantRecord, chrom: Chrom, region_searcher: SearchableSeqRegionsMap
):
    """Change `vcf_record` to be expressed relative to a different reference."""

    # Get index of personalised ref region containing the start of the `vcf_record`.
    region_index = region_searcher.bisect(chrom, vcf_record.pos, BisectTarget.PERS_REF)

    consumed_reference = 0  # Position in the inferred ref. sequence
    reference_length = len(vcf_record.ref)
    ref_seq_left = True

    # Build the rebased_ref as we traverse the regions,
    # using the vcf_record alt and pre-pend and/or post-pend to it if necessary.
    rebased_ref, rebased_alt = "", str(vcf_record.alts[0])

    # Let's rebase the position straight away
    first_region = region_searcher.get_region(chrom, region_index)

    # Case: hitting variant region. We rebase at the beginning of the variant region.
    if first_region.is_variant_region:
        rebased_pos = first_region.base_ref_start

        # We also straight away pre-pend any preceding variation relative to the base ref
        if vcf_record.pos > first_region.pers_ref_start:
            record_inset = vcf_record.pos - first_region.pers_ref_start
            rebased_alt = first_region.vcf_record_alt[:record_inset] + rebased_alt

    # Case: hitting non-variant region. We rebase at where the vcf_record starts, in base ref coordinates.
    else:
        rebased_pos = first_region.base_ref_start + (
            vcf_record.pos - first_region.pers_ref_start
        )

    while ref_seq_left:
        region = region_searcher.get_region(chrom, region_index)
        # Check how much of the vcf_record ref (inferred reference) can be consumed by the current region.
        # If the current region can consume at least what is left of the vcf_record ref, loop ends.
        # NOTE that region.length is 'overloaded': if a non-var region, it is the fixed interval between var regions.
        # If a var region, it is the inferred_vcf record's alt length (ref and alt lengths can differ).
        consumable = region.length - (
            vcf_record.pos + consumed_reference - region.pers_ref_start
        )

        if consumable >= (reference_length - consumed_reference):
            ref_seq_left = False
            to_consume = reference_length - consumed_reference
        else:
            to_consume = consumable

        if region.is_variant_region:
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
    if region.is_variant_region:
        cur_pos = vcf_record.pos + consumed_reference
        # The inset will be < 0 if there is a part of the (inferred vcf record's) alt which has not been
        inset = cur_pos - (region.pers_ref_start + region.length)
        if inset < 0:
            rebased_alt += region.vcf_record_alt[inset:]

    vcf_record = _modify_vcf_record(
        vcf_record, pos=rebased_pos, ref=rebased_ref, alts=[rebased_alt]
    )

    return vcf_record


def _add_contig_lines(disco_paths: DiscoverPaths):
    """
    If you don't add contig headers in new variant file before loading the records it contains,
    end up with downstream contig name bugs when writing rebased records.
    [TODO] I could not figure how to modify vcf file's headers only using pysam interface.
    """
    derived_vcf_contigs = VariantFile(disco_paths.discov_vcf_cortex).header.contigs

    if list(derived_vcf_contigs) == list():  # Need to add contig headers
        base_vcf_contigs = VariantFile(disco_paths.geno_vcf).header.contigs
        with open(disco_paths.discov_vcf_cortex) as f_in:
            all_lines = f_in.readlines()
        with open(disco_paths.discov_vcf_cortex, "w") as f_out:
            insert = True
            for line in all_lines:
                f_out.write(line)
                if (line[0] != "##") and insert:
                    for contig in base_vcf_contigs:
                        f_out.write(f"##contig=<ID={contig}>\n")
                    insert = False


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
    vcf_record: VariantRecord, inferred_sequence: str, var_unplaced_records: List[str]
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
