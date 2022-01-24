"""
Variant discovery against inferred personalised reference.
Works on a haploid personalised reference produced by gramtools `genotype`
Any newly discovered variants get rebased, ie expressed in the original prg
    reference coordinate space.
"""
from typing import List
import logging
import json
import sys

from pysam import VariantFile, VariantRecord

from gramtools.commands.paths import DiscoverPaths
from gramtools.commands.common import load_fasta
from gramtools.commands.genotype.seq_region_map import (
    Chrom,
    ChromSizes,
    SeqRegionsMap,
    SeqRegionMapper,
    SearchableSeqRegionsMap,
    BisectTarget,
)
from gramtools import py_cortex_api_message

log = logging.getLogger("gramtools")


def run(args):
    try:
        import cortex.calls as cortex
    except ModuleNotFoundError as err:
        print(
            "Cannot run discovery: missing cortex variant caller."
            f"{py_cortex_api_message}",
            file=sys.stderr,
        )
        exit(1)
    log.info("Start process: discover")
    disco_paths = DiscoverPaths(args.disco_dir, args.geno_dir, args.force)
    disco_paths.setup()

    enforce_genotyping_was_haploid(disco_paths)

    # call variants using `cortex`
    log.debug("Running cortex")
    cortex_args = {
        "reference_fasta": disco_paths.pers_ref,
        "reads_files": disco_paths.reads_files,
        "output_vcf_file_path": disco_paths.discov_vcf,
    }
    if hasattr(args, "mem_height"):  # Allows for low memory integration tests
        cortex_args["mem_height"] = getattr(args, "mem_height")
    cortex.run(**cortex_args)

    #  Convert coordinates in personalised reference space to coordinates in (original) prg space.
    log.debug("Rebasing vcf")
    rebased_vcf_records = _rebase_vcf(disco_paths)
    num_records = len(rebased_vcf_records)

    _dump_rebased_vcf(rebased_vcf_records, disco_paths)

    log.info(f"Found {num_records} variants. " f"Final vcf in {disco_paths.final_vcf}")
    log.info("End process: discover.")


def _rebase_vcf(disco_paths: DiscoverPaths, check_records=True):
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
    derived_records = VariantFile(disco_paths.discov_vcf).fetch()

    # Not loading genotype-produced rebasing map here, because it lacks the sequences
    chrom_sizes: ChromSizes = _load_contig_sizes_from_vcf(disco_paths.geno_vcf)
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
    template_vcf = VariantFile(disco_paths.discov_vcf)
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
    """
    Changes `vcf_record` to be expressed relative to a different reference.

    The algorithm is not trivial to understand- refer to the tests to understand
    expected inputs/outputs and to figures/text in the gramtools PhD thesis for details.

    Brief explanation:
       Notation:
          - base reference = reference on which to rebase
          - personalised reference = reference on which `vcf_record` variation lies
       Goal: get base reference sequence and position and new alt sequence for `vcf_record`
       Functioning:
          - We use a map (`region_searcher`) storing coordinates/sequences in both reference spaces
          - We initially bisect into the map at first position <= `vcf_record`'pos in
            personalised reference space. The algorithm then goes through the map until
            it reaches the end of `vcf_record`'s pos, constructing new sequence along
            the way.
          - If we start or end in a variant site in base reference space, we use all of
            the ref/alt sequences in the new ref/alt sequences, because we need to carry
            over variation that exists in the personalised reference
          - If we start or end in an invariant site, we only use sequence from/up to the
            `vcf_record`'s ref sequence to paste in to the new base reference sequence.
    """
    cur_region_index = region_searcher.bisect(
        chrom, vcf_record.pos, BisectTarget.PERS_REF
    )
    cur_region = region_searcher.get_region(chrom, cur_region_index)

    new_ref_seq = ""
    new_alt_seq = str(vcf_record.alts[0])
    cur_pers_ref_pos = vcf_record.pos

    pers_ref_end_pos = cur_pers_ref_pos + len(vcf_record.ref) - 1
    new_pos = cur_region.base_ref_start

    num_bases_past_first_region = cur_pers_ref_pos - cur_region.pers_ref_start
    if num_bases_past_first_region > 0:
        if cur_region.is_variant_region:
            new_alt_seq = (
                cur_region.vcf_record_alt[:num_bases_past_first_region] + new_alt_seq
            )
        else:
            new_pos += num_bases_past_first_region

    while cur_pers_ref_pos <= pers_ref_end_pos:
        cur_region = region_searcher.get_region(chrom, cur_region_index)
        cur_region_end = cur_region.pers_ref_start + cur_region.length - 1
        num_bases_past_last_region = max(cur_region_end - pers_ref_end_pos, 0)
        if cur_region.is_variant_region:
            new_ref_seq += cur_region.vcf_record_ref
        else:
            start_offset = cur_pers_ref_pos - vcf_record.pos
            end_offset = cur_region_end - vcf_record.pos - num_bases_past_last_region
            new_ref_seq += vcf_record.ref[start_offset : end_offset + 1]
        if num_bases_past_last_region > 0 and cur_region.is_variant_region:
            offset = cur_region.length - num_bases_past_last_region
            new_alt_seq = new_alt_seq + cur_region.vcf_record_alt[offset:]
        cur_pers_ref_pos = cur_region_end + 1
        cur_region_index += 1

    vcf_record = _modify_vcf_record(
        vcf_record, pos=new_pos, ref=new_ref_seq, alts=[new_alt_seq]
    )
    return vcf_record


def _load_contig_sizes_from_vcf(vcf_fname) -> ChromSizes:
    result = dict()
    contigs = VariantFile(vcf_fname).header.contigs
    for contig_ID, contig_properties in contigs.items():
        result[contig_ID] = contig_properties.length
    if len(result) == 0:
        raise ValueError(
            f"{vcf_fname} does not have 'contig' lines giving contig sizes"
        )
    return result


def _add_contig_lines(disco_paths: DiscoverPaths):
    """
    If you don't add contig headers in new variant file before loading the records it contains,
    end up with downstream contig name bugs when writing rebased records.
    [TODO] I could not figure how to modify vcf file's headers only using pysam interface.
    """
    derived_vcf_contigs = VariantFile(disco_paths.discov_vcf).header.contigs

    if list(derived_vcf_contigs) == list():  # Need to add contig headers
        base_vcf_contigs = VariantFile(disco_paths.geno_vcf).header.contigs
        with open(disco_paths.discov_vcf) as f_in:
            all_lines = f_in.readlines()
        with open(disco_paths.discov_vcf, "w") as f_out:
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
