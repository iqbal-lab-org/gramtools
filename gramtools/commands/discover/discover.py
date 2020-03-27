"""
Variant discovery against inferred personalised reference.
Works on a haploid personalised reference produced by gramtools `genotype`
Any newly discovered variants get rebased, ie expressed in the original prg
    reference coordinate space.
"""
from typing import List, Dict
import bisect
import logging
import collections
import json

from pysam import VariantFile, VariantRecord
import cortex
from Bio import SeqIO
from pathlib import Path

from gramtools.commands.paths import DiscoverPaths

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

    chrom_sizes = list()
    with open(disco_paths.pers_ref) as genome_file:
        for seq_record in SeqIO.parse(genome_file, "fasta"):
            chrom_sizes.append(len(seq_record.seq))

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


## Rebase a vcf so that it uses same reference as base_vcf.
#
# Here is the use case. We have:
#  discovery.vcf                   personalised_ref.vcf
#   |                               |
#  personalised_ref.fasta          base_ref.fasta
#
# And we want:
#  discovery.vcf
#   |
#  base_ref.fasta
#
# In the context of gramtools:
#  * `discovery.vcf` is the vcf produced by cortex, ie variant calls against personalised reference.
# * `personalised_ref` is the personalised reference; it is the original prg with inferred variant alleles instead.
# * `personalised_ref.vcf` is the vcf describing the personalised reference with respect to the original prg.
#  * `base_ref` is the original prg. The first allele of each variant site in the prg is taken for making this.
#
# Thus what rebasing achieves is to express the variant calls in terms of the original prg coordinates.
def _rebase_vcf(disco_paths: DiscoverPaths, chrom_sizes, check_records=True):
    base_vcf_file_path = disco_paths.geno_vcf
    var_unplaced_records = []

    if check_records:
        inferred_refs = load_multi_fasta(disco_paths.pers_ref)

    base_records = VariantFile(base_vcf_file_path)
    derived_records = VariantFile(disco_paths.discov_vcf_cortex)

    # Given a variant call position w.r.t the personalised reference, we will need to know how this position maps to the original prg.
    flagged_personalised_regions = _flag_personalised_reference_regions(
        base_records, chrom_sizes
    )

    new_vcf_records = []
    for vcf_record in derived_records.fetch():
        chrom_key = vcf_record.chrom

        # Check the variant call is properly placed by the variant caller
        if check_records:
            position = vcf_record.pos
            inferred_sequence = inferred_refs[chrom_key]
            if (
                len(inferred_sequence) < position
                or vcf_record.ref
                != inferred_sequence[position - 1 : position - 1 + len(vcf_record.ref)]
            ):
                var_unplaced_records.append(str(vcf_record))
                continue

        regions_list = flagged_personalised_regions.get(chrom_key, None)
        assert (
            regions_list is not None
        ), "Error: ref ID {} in vcf record {} is not present in vcf from infer {}".format(
            chrom_key, vcf_record, base_vcf_file_path
        )

        new_record = _rebase_vcf_record(vcf_record, regions_list)

        # If var caller calls reverses `genotype`'s call, ignore
        if not all([new_record.ref == alt for alt in new_record.alts]):
            new_vcf_records.append(new_record)

    if len(var_unplaced_records) > 0:
        log.warning(
            "The following new variant records were skipped: {} \n"
            "Reasons: POS is not a valid position in the inferred reference sequence,"
            "or REF does not coincide with inferred reference sequence.".format(
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


class _Region:
    """Mapping between vcf records in two coordinate spaces
    """

    def __init__(
        self, base_POS, inf_POS, length, vcf_record_REF=None, vcf_record_ALT=None
    ):
        # Start coordinates
        self.base_POS = base_POS
        self.inf_POS = inf_POS

        self.length = length
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


_Regions = List[_Region]
_Regions_Map = Dict[str, _Regions]


def _flag_personalised_reference_regions(
    base_records: VariantFile, chrom_sizes
) -> _Regions_Map:
    """
    There are two types of `Regions`: those inside variant sites in the prg, and those outside.
    Important requirements for vcf:
        * the records are sorted by POS for a given ref ID, and all records referring to the same ref ID are contiguous.
        * the chrom_sizes are ordered in same order as the ref IDs in the records.
    The first point is explicitly checked using assertions.
    """
    all_personalised_ref_regions = {}  # Maps each CHROM to a list of _Regions
    chrom_positions = {}  # Maps each CHROM to a base_ref and an inferred_ref position
    num_chroms = -1
    prev_chrom_key = None
    prev_record = None

    for vcf_record in base_records.fetch():
        chrom_key = vcf_record.chrom
        if chrom_key not in all_personalised_ref_regions:
            all_personalised_ref_regions[chrom_key] = []
            chrom_positions[chrom_key] = [1, 1]

            # We enter a new ref ID; let's make sure we capture the end of the previous ref ID
            if num_chroms >= 0:
                base_ref_pos, inf_ref_pos = (
                    chrom_positions[prev_chrom_key][0],
                    chrom_positions[prev_chrom_key][1],
                )
                if base_ref_pos <= chrom_size:
                    regions_list = all_personalised_ref_regions[prev_chrom_key]
                    non_var_region = _Region(
                        base_POS=base_ref_pos,
                        inf_POS=inf_ref_pos,
                        length=chrom_size - base_ref_pos + 1,
                    )
                    regions_list.append(non_var_region)

            num_chroms += 1
            chrom_size = chrom_sizes[num_chroms]

        else:
            # Enforce ref ID contiguity
            assert (
                chrom_key == prev_chrom_key
            ), "Ref IDs not contiguous: {} and {} interspersed".format(
                chrom_key, prev_chrom_key
            )
            # Enforce position sortedness
            assert (
                vcf_record.pos > prev_record.pos
            ), "Records not in increasing POS order: {} and {}".format(
                prev_record, vcf_record
            )
        regions_list = all_personalised_ref_regions[chrom_key]
        base_ref_pos, inf_ref_pos = (
            chrom_positions[chrom_key][0],
            chrom_positions[chrom_key][1],
        )

        # Build non-variant region
        if vcf_record.pos > base_ref_pos:
            region_length = vcf_record.pos - base_ref_pos
            if len(regions_list) == 0 or regions_list[-1].is_site:
                non_var_region = _Region(
                    base_POS=base_ref_pos, inf_POS=inf_ref_pos, length=region_length
                )
                regions_list.append(non_var_region)
            else:  # We saw a site, but it was not used; so, we only need to merge this region with prev region.
                regions_list[-1].length += region_length
            base_ref_pos += region_length
            inf_ref_pos += region_length

        # Build variant region
        picked_alleles = vcf_record.samples[0]["GT"]

        if set(picked_alleles) == {None}:
            picked_allele = 0  # null GT, take REF allele
        else:
            picked_allele = picked_alleles[0]

        # Only make a var region if the inference procedure did not pick the REF.
        if picked_allele != 0:
            var_region = _Region(
                base_POS=base_ref_pos,
                inf_POS=inf_ref_pos,
                length=len(vcf_record.ALT[picked_allele - 1]),
                vcf_record_REF=vcf_record.REF,
                vcf_record_ALT=str(vcf_record.ALT[picked_allele - 1]),
            )

            regions_list.append(var_region)
            base_ref_pos += len(vcf_record.REF)
            inf_ref_pos += var_region.length

        chrom_positions[chrom_key] = [base_ref_pos, inf_ref_pos]
        prev_chrom_key = chrom_key
        prev_record = vcf_record

    # End game: deal with the last non-var region if there is one.
    if len(all_personalised_ref_regions) == 0:
        raise ValueError("No records in provided vcf.")

    if base_ref_pos <= chrom_size:
        non_var_region = _Region(
            base_POS=base_ref_pos,
            inf_POS=inf_ref_pos,
            length=chrom_size - base_ref_pos + 1,
        )
        regions_list.append(non_var_region)

    return all_personalised_ref_regions


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

    # we will build the rebased_REF as we traverse the regions. We will use the vcf_record ALT and pre-pend and/or post-pend to it if necessary.
    rebased_REF, rebased_ALT = "", str(vcf_record.alts[0])

    # Let's rebase the position straight away
    first_region = personalised_reference_regions[region_index]

    # Case: hitting variant region. We rebase at the beginning of the variant region.
    if first_region.is_site:
        rebased_POS = first_region.base_POS

        # We also straight away pre-pend any preceding variation relative to the base REF
        if vcf_record.pos > first_region.inf_POS:
            record_inset = vcf_record.pos - first_region.inf_POS
            rebased_ALT = first_region.vcf_record_ALT[:record_inset] + rebased_ALT

    # Case: hitting non-variant region. We rebase at where the vcf_record starts, in base ref coordinates.
    else:
        rebased_POS = first_region.base_POS + (vcf_record.pos - first_region.inf_POS)

    while ref_seq_left:
        region = personalised_reference_regions[region_index]
        # Check how much of the vcf_record REF (inferred reference) can be consumed by the current region.
        # If the current region can consume at least what is left of the vcf_record REF, loop ends.
        # NOTE that region.length is 'overloaded': if a non-var region, it is the fixed interval between var regions.
        # If a var region, it is the inferred_vcf record's ALT length (REF and ALT lengths can differ).
        consumable = region.length - (
            vcf_record.pos + consumed_reference - region.inf_POS
        )

        if consumable >= (reference_length - consumed_reference):
            ref_seq_left = False
            to_consume = reference_length - consumed_reference
        else:
            to_consume = consumable

        if region.is_site:
            rebased_REF += region.vcf_record_REF

        else:
            # We can use the vcf_record's REF, as that is also the base REF sequence- because we are in a non-variant site.
            rebased_REF += vcf_record.ref[
                consumed_reference : consumed_reference + to_consume
            ]

        consumed_reference += to_consume
        region_index += 1

    assert consumed_reference == len(vcf_record.ref)

    # Deal with the last region: post-pend any sequence in ALT record if we finish in a variant site.
    if region.is_site:
        cur_pos = vcf_record.pos + consumed_reference
        # The inset will be < 0 if there is a part of the (inferred vcf record's) ALT which has not been
        inset = cur_pos - (region.inf_POS + region.length)
        if inset < 0:
            rebased_ALT += region.vcf_record_ALT[inset:]

    vcf_record = _modify_vcf_record(
        vcf_record, POS=rebased_POS, REF=rebased_REF, ALT=[rebased_ALT]
    )

    return vcf_record


## Return the marked `_Region` containing the position referred to by `vcf_record`.
def _find_start_region_index(vcf_record: VariantRecord, marked_regions):
    record_start_POS = vcf_record.pos

    # Binary search on a sorted array. The '<' comparator is overloaded in `_Region` definition to allow the search.
    # The produced index is:
    # * If a region with start POS == record_start_POS is found, returns that index (specifically, the first occurrence)
    # * If such a region is not found, returns the index of the region with start POS **strictly** greater than record_start_POS
    # IN THAT CASE, we need to return the region index one smaller; so that the region's interval includes record_start_POS.
    region_index = bisect.bisect_left(marked_regions, record_start_POS)

    #  Case: `record_start_index` is larger than any _Region start.
    if region_index > len(marked_regions) - 1:
        region_index = len(marked_regions) - 1

    selected_region = marked_regions[region_index]

    # Case: record_start_POS was not found exactly.
    if selected_region.inf_POS > record_start_POS:
        region_index -= 1

    return region_index


## Produce a list of vcf records parsed from `vcf_file_path`.
def _load_records(vcf_file_path: Path):
    vcf_reader = vcf.Reader(filename=str(vcf_file_path))
    all_vcf_records = list(vcf_reader)
    return all_vcf_records


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


## Generator to always yield the first allele (REF) of a variant site.
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
