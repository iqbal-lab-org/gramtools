from typing import List, Dict, Iterable, Union, Callable
from enum import Enum, auto
import json
from pathlib import Path

from pysam import VariantRecord

VariantRecords = Iterable[VariantRecord]
Chrom = str
ChromSizes = Dict[Chrom, int]


class SeqRegion:
    """Mapping between vcf records in two coordinate spaces
    """

    def __init__(
        self,
        base_ref_start: int,
        pers_ref_start: int,
        length: int,
        vcf_record_ref: Union[str, None] = None,
        vcf_record_alt: Union[str, None] = None,
    ):
        # Start coordinates
        self.base_ref_start = base_ref_start
        self.pers_ref_start = pers_ref_start

        self.vcf_record_ref = vcf_record_ref
        self.vcf_record_alt = vcf_record_alt

        if self.vcf_record_alt is not None and length is not None:
            if length != len(self.vcf_record_alt):
                raise ValueError(
                    f"{length} must be length of {vcf_record_alt} when both are provided."
                )
        self.length = length

    @property
    def is_variant_region(self):
        return self.vcf_record_ref is not None

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __repr__(self):
        return str(self.__dict__)

    def to_json(self, dump_sequences: bool = True) -> Dict:
        dumped_attrs = self.__dict__.copy()
        if not dump_sequences:
            dumped_attrs["vcf_record_ref"] = dumped_attrs["vcf_record_alt"] = None
        for attr in self.__dict__:
            if dumped_attrs[attr] is None:
                dumped_attrs.pop(attr)
        return {"SeqRegion": dumped_attrs}

    @staticmethod
    def from_json(dct: Dict) -> Union["SeqRegion", Dict]:
        if "SeqRegion" in dct:
            return SeqRegion(**dct["SeqRegion"])
        return dct


SeqRegions = List[SeqRegion]
SeqRegionsMap = Dict[Chrom, SeqRegions]


class _PosTracker:
    def __init__(self, base_ref_pos: int, pers_ref_pos: int):
        self.base_ref_pos = base_ref_pos
        self.pers_ref_pos = pers_ref_pos


class SeqRegionMapper:
    def __init__(self, base_records: VariantRecords, chrom_sizes: ChromSizes):
        self.chrom_sizes = chrom_sizes
        self.map: SeqRegionsMap = dict()
        self.pos_trackers: Dict[Chrom, _PosTracker] = {}

        prev_chrom_key, prev_record = None, None

        for record in base_records:
            chrom_key = record.chrom
            if chrom_key not in self.map:
                self._new_chrom(chrom_key, prev_chrom_key)
            else:
                self._enforce_contiguity(chrom_key, prev_chrom_key, record, prev_record)

            base_pos = self.pos_trackers[chrom_key].base_ref_pos
            if record.pos > base_pos:  # Build non-variant region
                region_length = record.pos - base_pos
                self._add_invariant_region(chrom_key, region_length)

            self._add_variant_region(chrom_key, record)
            prev_chrom_key = chrom_key
            prev_record = record

        if len(self.map) == 0:
            raise ValueError("No records in provided vcf.")

        chrom_size = self.chrom_sizes[chrom_key]
        base_pos = self.pos_trackers[chrom_key].base_ref_pos
        if base_pos <= chrom_size:
            self._add_invariant_region(chrom_key, chrom_size - base_pos + 1)

        self.map_invariant_chroms()

    def get_map(self):
        return self.map

    @property
    def num_processed_chroms(self):
        return len(self.map.keys())

    def map_invariant_chroms(self):
        """
        The PRG may contain contigs with no variation.
        We need to map those so that variants found against them are recognised
        """
        for chrom in self.chrom_sizes:
            if chrom not in self.map:
                self.map[chrom] = [SeqRegion(1, 1, self.chrom_sizes[chrom])]

    def _new_chrom(self, chrom_key, prev_chrom_key):
        # New chrom; make sure capture the end of the previous chrom
        if self.num_processed_chroms > 0:
            prev_base_pos = self.pos_trackers[prev_chrom_key].base_ref_pos
            prev_chrom_size = self.chrom_sizes[prev_chrom_key]
            if prev_base_pos <= prev_chrom_size:
                region_length = prev_chrom_size - prev_base_pos + 1
                self._add_invariant_region(prev_chrom_key, region_length)

        self.map[chrom_key] = list()
        self.pos_trackers[chrom_key] = _PosTracker(1, 1)

    def _add_invariant_region(self, chrom_key, region_length: int):
        ref_positions = self.pos_trackers[chrom_key]

        focal_regions: SeqRegions = self.map[chrom_key]
        perform_extension = (
            len(focal_regions) > 0 and not focal_regions[-1].is_variant_region
        )
        if perform_extension:
            # Case: REF called variant region follows, or is followed by, invariant region
            focal_regions[-1].length += region_length
        else:
            non_var_region = SeqRegion(
                base_ref_start=ref_positions.base_ref_pos,
                pers_ref_start=ref_positions.pers_ref_pos,
                length=region_length,
            )
            self.map[chrom_key].append(non_var_region)
        ref_positions.base_ref_pos += region_length
        ref_positions.pers_ref_pos += region_length

    def _add_variant_region(self, chrom_key, vcf_record: VariantRecord):
        ref_positions = self.pos_trackers[chrom_key]
        picked_alleles = vcf_record.samples[0]["GT"]
        if set(picked_alleles) == {None}:
            picked_allele = 0  # null GT, take ref allele
        else:
            picked_allele = picked_alleles[0]

        if picked_allele != 0:
            var_region = SeqRegion(
                base_ref_start=ref_positions.base_ref_pos,
                pers_ref_start=ref_positions.pers_ref_pos,
                length=len(vcf_record.alts[picked_allele - 1]),
                vcf_record_ref=vcf_record.ref,
                vcf_record_alt=str(vcf_record.alts[picked_allele - 1]),
            )

            self.map[chrom_key].append(var_region)
            ref_positions.base_ref_pos += len(vcf_record.ref)
            ref_positions.pers_ref_pos += var_region.length
        else:
            self._add_invariant_region(chrom_key, len(vcf_record.ref))

    def _enforce_contiguity(self, chrom, prev_chrom, vcf_rec, prev_vcf_rec):
        # Enforce ref ID contiguity and position sortedness
        assert (
            chrom == prev_chrom
        ), f"Ref IDs not contiguous: {chrom} and {prev_chrom} interspersed"
        assert (
            vcf_rec.pos > prev_vcf_rec.pos
        ), f"Records not in increasing pos order: {prev_vcf_rec} and {vcf_rec}"


class BisectTarget(Enum):
    BASE_REF = auto()
    PERS_REF = auto()


class SearchableSeqRegionsMap:
    def __init__(self, map: SeqRegionsMap):
        self._map = map

    def bisect(self, chrom: Chrom, pos: int, mode: BisectTarget):
        if mode not in BisectTarget:
            raise ValueError(f"mode argument should be of type {BisectTarget}")
        regions: SeqRegions = self._map[chrom]
        if mode is BisectTarget.BASE_REF:
            return self._bisect(regions, pos, self.base_ref_gt)
        else:
            return self._bisect(regions, pos, self.pers_ref_gt)

    def get_region(self, chrom: Chrom, region_index: int):
        return self._map[chrom][region_index]

    @staticmethod
    def base_ref_gt(region: SeqRegion, pos: int):
        return region.base_ref_start > pos

    @staticmethod
    def pers_ref_gt(region: SeqRegion, pos: int):
        return region.pers_ref_start > pos

    def _bisect(
        self,
        regions: SeqRegions,
        pos: int,
        cmp_function: Callable[[SeqRegion, int], bool],
    ):
        """
        First, does right bisection: find the index of the element in `regions`
            for which all preceding elements have start position <= `pos`,
            and all other elements have start position > `pos`
        Then returns that - 1 so that we get the region which has start position <= `pos`
        """
        lo = 0
        hi = len(regions)
        while lo < hi:
            mid = (lo + hi) // 2
            if cmp_function(regions[mid], pos):
                hi = mid
            else:
                lo = mid + 1
        return lo - 1

    def __eq__(self, other: "SearchableSeqRegionsMap"):
        return self._map == other._map

    class JsonEncode(json.JSONEncoder):
        dump_sequences: bool

        def default(self, obj):
            if isinstance(obj, SeqRegion):
                return obj.to_json(self.dump_sequences)
            return json.JSONEncoder.default(self, obj)

    def dump_to(self, fname: Path, dump_sequences=True) -> None:
        self.JsonEncode.dump_sequences = dump_sequences
        with fname.open("w") as fout:
            json.dump(self._map, fout, cls=self.JsonEncode)

    @staticmethod
    def load_from(fname: Path) -> "SearchableSeqRegionsMap":
        with fname.open("r") as fin:
            loaded_map = json.load(fin, object_hook=SeqRegion.from_json)
            return SearchableSeqRegionsMap(loaded_map)
