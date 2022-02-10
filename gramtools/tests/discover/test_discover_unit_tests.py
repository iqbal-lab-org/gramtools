from typing import List
from unittest import TestCase, main as unittest_main

from gramtools.tests.mocks import _MockVcfRecord
from gramtools.commands.discover import discover
from gramtools.commands.genotype.seq_region_map import (
    SeqRegionMapper,
    SearchableSeqRegionsMap,
)

MockVcfRecords = List[_MockVcfRecord]


def make_map(
    base_records: MockVcfRecords, chrom_sizes: List[int]
) -> SearchableSeqRegionsMap:
    names = [f"chr{i}" for i in range(len(chrom_sizes))]
    named_chroms = dict(zip(names, chrom_sizes))
    region_map = SeqRegionMapper(base_records, named_chroms).get_map()
    return SearchableSeqRegionsMap(region_map)


def run_rebase(
    discov_record: _MockVcfRecord, base_records: MockVcfRecords, chrom_sizes: List[int]
) -> _MockVcfRecord:
    region_searcher = make_map(base_records, chrom_sizes)
    new_vcf_record = discover._rebase_vcf_record(
        discov_record, discov_record.chrom, region_searcher
    )
    return new_vcf_record


class TestRebasing_InvalidRebasing(TestCase):
    def test_rebasing_in_unknown_chromosome_fails(self):
        base_records = [_MockVcfRecord(pos=2, ref="T", alts=["G"], chrom="chr0")]
        discov_record = _MockVcfRecord(pos=3, ref="C", alts=["G"], chrom="chr1")
        with self.assertRaises(KeyError):
            result = run_rebase(discov_record, base_records, chrom_sizes=[5])

    def test_rebasing_with_too_short_chromosome_fails(self):
        # base sequence:      AA ATCTA
        # derived sequence:   T  ATCTA
        base_records = [_MockVcfRecord(pos=1, ref="AA", alts=["T"], chrom="chr0")]
        discov_record = _MockVcfRecord(pos=4, ref="C", alts=["G"], chrom="chr0")

        with self.assertRaises(IndexError):
            result = run_rebase(discov_record, base_records, chrom_sizes=[3])


class TestRebasing_ContainedRecords(TestCase):
    """
    Contained means the vcf record intersects with a single
    region in the rebasing map
    """

    def test_VariantInInvariantChromosome(self):
        base_records = [_MockVcfRecord(pos=2, ref="TAT", alts=["G"], chrom="chr0")]
        discov_record = _MockVcfRecord(pos=1, ref="TT", alts=["GA"], chrom="chr1")
        result = run_rebase(discov_record, base_records, chrom_sizes=[7, 5])

        self.assertEqual(discov_record, result)

    def test_VariantCoveringAllOfInvariantRegion(self):
        # base sequence:      AA ATATA
        # derived sequence:   T  ATATA
        base_records = [_MockVcfRecord(pos=1, ref="AA", alts=["T"], chrom="chr0")]
        discov_record = _MockVcfRecord(pos=2, ref="ATATA", alts=["C"], chrom="chr0")

        expected = _MockVcfRecord(pos=3, ref="ATATA", alts=["C"], chrom="chr0")
        result = run_rebase(discov_record, base_records, chrom_sizes=[7])
        self.assertEqual(expected, result)

    def test_VariantCoveringPartOfInvariantRegion(self):
        # base sequence:      AA ATCTA
        # derived sequence:   T  ATCTA
        base_records = [_MockVcfRecord(pos=1, ref="AA", alts=["T"], chrom="chr0")]
        discov_record = _MockVcfRecord(pos=4, ref="C", alts=["G"], chrom="chr0")

        expected = _MockVcfRecord(pos=5, ref="C", alts=["G"], chrom="chr0")
        result = run_rebase(discov_record, base_records, chrom_sizes=[7])
        self.assertEqual(expected, result)

    def test_VariantCoveringAllOfVariantRegion(self):
        # base sequence:      - TAT ---
        # derived sequence:   - G   ---
        base_records = [_MockVcfRecord(pos=2, ref="TAT", alts=["G"], chrom="chr0")]
        discov_record = _MockVcfRecord(pos=3, ref="G", alts=["C"], chrom="chr0")

        result = run_rebase(discov_record, base_records, chrom_sizes=[7])
        expected = _MockVcfRecord(pos=5, ref="G", alts=["C"], chrom="chr0")
        self.assertEqual(expected, result)

    def test_VariantCoveringPartOfVariantRegion(self):
        # base sequence:      - TAAAT ---
        # derived sequence:   - TAT   ---
        base_records = [_MockVcfRecord(pos=2, ref="TAAAT", alts=["TAT"], chrom="chr0")]
        discov_record = _MockVcfRecord(pos=3, ref="A", alts=["C"], chrom="chr0")

        result = run_rebase(discov_record, base_records, chrom_sizes=[9])
        expected = _MockVcfRecord(pos=2, ref="TAAAT", alts=["TCT"], chrom="chr0")
        self.assertEqual(expected, result)


class TestRebasing_OverlappingRecords(TestCase):
    """
    Overlapping means the vcf record spans >1 regions in the rebasing map
    """

    def test_VariantOverlapsTwoRegions_AllOfVarPartOfInvar(self):
        # base sequence:      AAA  AGA A
        # derived sequence:   TTTT AGA C
        # overlap:            ---- --
        base_records = [
            _MockVcfRecord(pos=1, ref="AAA", alts=["TTTT"], chrom="chr0"),
            _MockVcfRecord(pos=7, ref="A", alts=["C"], chrom="chr0"),
        ]
        discov_record = _MockVcfRecord(
            pos=1, ref="TTTTAG", alts=["TATTAC"], chrom="chr0"
        )
        result = run_rebase(discov_record, base_records, chrom_sizes=[7])
        expected = _MockVcfRecord(pos=1, ref="AAAAG", alts=["TATTAC"], chrom="chr0")

        self.assertEqual(expected, result)

    def test_VariantOverlapsTwoRegions_PartOfVarAllOfInvar(self):
        # base sequence:      AAA  AGA A
        # derived sequence:   TTTT AGA C
        # overlap:              -- ---
        base_records = [
            _MockVcfRecord(pos=1, ref="AAA", alts=["TTTT"], chrom="chr0"),
            _MockVcfRecord(pos=7, ref="A", alts=["C"], chrom="chr0"),
        ]
        discov_record = _MockVcfRecord(pos=3, ref="TTAGA", alts=["TATGA"], chrom="chr0")
        result = run_rebase(discov_record, base_records, chrom_sizes=[7])
        expected = _MockVcfRecord(pos=1, ref="AAAAGA", alts=["TTTATGA"], chrom="chr0")

        self.assertEqual(expected, result)

    def test_VariantOverlapsTwoRegions_PartOfInvarAllOfVar(self):
        # base sequence:      AAA  AGA A
        # derived sequence:   TTTT AGA C
        # overlap:                  -- -
        base_records = [
            _MockVcfRecord(pos=1, ref="AAA", alts=["TTTT"], chrom="chr0"),
            _MockVcfRecord(pos=7, ref="A", alts=["C"], chrom="chr0"),
        ]
        discov_record = _MockVcfRecord(pos=6, ref="GAC", alts=["AAT"], chrom="chr0")
        result = run_rebase(discov_record, base_records, chrom_sizes=[7])
        expected = _MockVcfRecord(pos=5, ref="GAA", alts=["AAT"], chrom="chr0")

        self.assertEqual(expected, result)

    def test_VariantOverlapsThreeRegions_VarThenInvarThenVar_FullSpan(self):
        # base sequ: T TAT GGG T     ATTTT
        # secondary: T GG  GGG TCTGT ATTTT
        # overlap:     --  --- -----
        base_records = [
            _MockVcfRecord(pos=2, ref="TAT", alts=["GG"], chrom="chr0"),
            _MockVcfRecord(pos=8, ref="T", alts=["TCTGT"], chrom="chr0"),
        ]
        discov_record = _MockVcfRecord(
            pos=2, ref="GGGGGTCTGT", alts=["GAGAGTCAGT"], chrom="chr0"
        )
        result = run_rebase(discov_record, base_records, chrom_sizes=[13])
        expected = _MockVcfRecord(
            pos=2, ref="TATGGGT", alts=["GAGAGTCAGT"], chrom="chr0"
        )
        self.assertEqual(expected, result)

    def test_VariantOverlapsThreeRegions_VarThenInvarThenVar_PartialSpan(self):
        # base sequ: T TAT GGG T     ATTTT
        # secondary: T GG  GGG TCTGT ATTTT
        # overlap:      -  --- ---
        base_records = [
            _MockVcfRecord(pos=2, ref="TAT", alts=["GG"], chrom="chr0"),
            _MockVcfRecord(pos=8, ref="T", alts=["TCTGT"], chrom="chr0"),
        ]

        discov_record = _MockVcfRecord(
            pos=3, ref="GGGGTCT", alts=["ACCCTCA"], chrom="chr0"
        )
        result = run_rebase(discov_record, base_records, chrom_sizes=[13])
        expected = _MockVcfRecord(
            pos=2, ref="TATGGGT", alts=["GACCCTCAGT"], chrom="chr0"
        )
        self.assertEqual(expected, result)

    def test_VariantOverlapsThreeRegions_InVarThenVarThenInvar_FullSpan(self):
        # base sequ: T TAT GGG T     ATTTT
        # secondary: T GG  GGG TCTGT ATTTT
        # overlap:         --- ----- -----

        base_records = [
            _MockVcfRecord(pos=2, ref="TAT", alts=["GG"], chrom="chr0"),
            _MockVcfRecord(pos=8, ref="T", alts=["TCTGT"], chrom="chr0"),
        ]
        discov_record = _MockVcfRecord(
            pos=4, ref="GGGTCTGTATTTT", alts=["GCGTCAGTATTCT"], chrom="chr0"
        )
        result = run_rebase(discov_record, base_records, chrom_sizes=[13])
        expected = _MockVcfRecord(
            pos=5, ref="GGGTATTTT", alts=["GCGTCAGTATTCT"], chrom="chr0"
        )
        self.assertEqual(expected, result)

    def test_VariantOverlapsThreeRegions_InVarThenVarThenInvar_PartialSpan(self):
        # base sequ: T TAT GGG T     ATTTT
        # secondary: T GG  GGG TCTGT ATTTT
        # overlap:          -- ----- --

        base_records = [
            _MockVcfRecord(pos=2, ref="TAT", alts=["GG"], chrom="chr0"),
            _MockVcfRecord(pos=8, ref="T", alts=["TCTGT"], chrom="chr0"),
        ]
        discov_record = _MockVcfRecord(pos=5, ref="GGTCTGTAT", alts=["T"], chrom="chr0")
        result = run_rebase(discov_record, base_records, chrom_sizes=[13])
        expected = _MockVcfRecord(pos=6, ref="GGTAT", alts=["T"], chrom="chr0")
        self.assertEqual(expected, result)


if __name__ == "__main__":
    unittest_main()
