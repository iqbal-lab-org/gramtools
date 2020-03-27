from typing import List, Dict, Iterable
from pysam import VariantFile, VariantRecord

VariantRecords = Iterable[VariantRecord]


class _Region:
    """Mapping between vcf records in two coordinate spaces
    """

    def __init__(
        self, base_pos, inf_pos, length, vcf_record_ref=None, vcf_record_alt=None
    ):
        # Start coordinates
        self.base_pos = base_pos
        self.inf_pos = inf_pos

        self.length = length
        self.vcf_record_ref = vcf_record_ref
        self.vcf_record_alt = vcf_record_alt

    @property
    def is_site(self):
        return self.vcf_record_ref is not None

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __repr__(self):
        return str(self.__dict__)

    def __lt__(self, other):
        if isinstance(other, _Region):
            return self.inf_pos < other.inf_pos
        else:
            return self.inf_pos < other


_Regions = List[_Region]
_Regions_Map = Dict[str, _Regions]


class _ref_positions:
    def __init__(self, base_pos, derived_pos):
        self.base_pos = base_pos
        self.derived_pos = derived_pos


class RegionMapper:
    def __init__(self, base_records: VariantRecords, chrom_sizes):
        self.base_recs = base_records
        self.chrom_sizes = chrom_sizes
        self.all_regions = {}  # Maps each CHROM to a list of _Regions
        self.chrom_pos: Dict[str, _ref_positions] = {}

    def get_mapped(self):
        if len(self.all_regions) > 0:
            return self.all_regions

        prev_chrom_key, prev_record = None, None

        for record in self.base_recs:
            chrom_key = record.chrom
            # Switch contigs
            if chrom_key not in self.all_regions:
                self._new_chrom(chrom_key, prev_chrom_key)
            else:
                self._enforce_contiguity(chrom_key, prev_chrom_key, record, prev_record)

            # Build non-variant region
            base_pos = self.chrom_pos[chrom_key].base_pos
            if record.pos > base_pos:
                region_length = record.pos - base_pos
                assert (
                    len(self.all_regions[chrom_key]) == 0
                    or self.all_regions[chrom_key][-1].is_site
                )
                self._add_invariant_region(chrom_key, region_length)

            self._add_variant_region(chrom_key, record)
            prev_chrom_key = chrom_key
            prev_record = record

        if len(self.all_regions) == 0:
            raise ValueError("No records in provided vcf.")

        base_pos = self.chrom_pos[chrom_key].base_pos
        if base_pos <= self.chrom_size:
            self._add_invariant_region(chrom_key, self.chrom_size - base_pos + 1)

        return self.all_regions

    @property
    def num_chroms(self):
        return len(self.all_regions.keys())

    @property
    def chrom_size(self):
        return self.chrom_sizes[self.num_chroms - 1]

    def _new_chrom(self, chrom_key, prev_chrom_key):
        # New chrom; make sure capture the end of the previous chrom
        if self.num_chroms > 0:
            prev_base_pos = self.chrom_pos[prev_chrom_key].base_pos
            if prev_base_pos <= self.chrom_size:
                region_length = self.chrom_size - prev_base_pos + 1
                self._add_invariant_region(prev_chrom_key, region_length)

        self.all_regions[chrom_key] = []
        self.chrom_pos[chrom_key] = _ref_positions(1, 1)

    def _add_invariant_region(self, chrom_key, region_length: int):
        ref_positions = self.chrom_pos[chrom_key]
        non_var_region = _Region(
            base_pos=ref_positions.base_pos,
            inf_pos=ref_positions.derived_pos,
            length=region_length,
        )
        self.all_regions[chrom_key].append(non_var_region)
        ref_positions.base_pos += region_length
        ref_positions.derived_pos += region_length

    def _add_variant_region(self, chrom_key, vcf_record: VariantRecord):
        ref_positions = self.chrom_pos[chrom_key]
        picked_alleles = vcf_record.samples[0]["GT"]
        if set(picked_alleles) == {None}:
            picked_allele = 0  # null GT, take ref allele
        else:
            picked_allele = picked_alleles[0]

        if picked_allele != 0:
            var_region = _Region(
                base_pos=ref_positions.base_pos,
                inf_pos=ref_positions.derived_pos,
                length=len(vcf_record.alts[picked_allele - 1]),
                vcf_record_ref=vcf_record.ref,
                vcf_record_alt=str(vcf_record.alts[picked_allele - 1]),
            )

            self.all_regions[chrom_key].append(var_region)
            ref_positions.derived_pos += var_region.length
        else:
            ref_positions.derived_pos += len(vcf_record.ref)

        ref_positions.base_pos += len(vcf_record.ref)

    def _enforce_contiguity(self, chrom, prev_chrom, vcf_rec, prev_vcf_rec):
        # Enforce ref ID contiguity and position sortedness
        assert (
            chrom == prev_chrom
        ), f"Ref IDs not contiguous: {chrom} and {prev_chrom} interspersed"
        assert (
            vcf_rec.pos > prev_vcf_rec.pos
        ), f"Records not in increasing pos order: {prev_vcf_rec} and {vcf_rec}"
