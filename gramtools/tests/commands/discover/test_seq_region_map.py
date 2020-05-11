from unittest import TestCase

from gramtools.tests.mocks import _MockVcfRecord
from gramtools.commands.discover.seq_region_map import SeqRegion
from gramtools.commands.discover.discover import RegionMapper, _find_start_region_index


class TestRegionMapping(TestCase):
    """
    Note: as stated in the original code, the length of a region is the length of the alt string if we have a vcf record.
    """

    def test_SingleBaseAlt_CorrectRegion(self):
        # base sequence:      T TAT CGG
        # derived sequence:   T G   CGG
        base_records = [_MockVcfRecord(pos=2, ref="TAT", alts=["G"])]

        chrom_sizes = {"JAC": 7}

        result = RegionMapper(base_records, chrom_sizes).get_map()

        expected = [
            SeqRegion(base_pos=1, inf_pos=1, length=1),
            SeqRegion(
                base_pos=2,
                inf_pos=2,
                length=1,
                vcf_record_ref="TAT",
                vcf_record_alt="G",
            ),
            SeqRegion(base_pos=5, inf_pos=3, length=3),
        ]
        self.assertEqual(expected, result["JAC"])

    def test_ref_call_produces_invariant_region_only(self):
        # base sequence:      T TAT CGG
        # derived sequence:   ^^^^^^^^^
        base_records = [
            _MockVcfRecord(pos=2, ref="TAT", alts=["G"], samples=[{"GT": [0]}])
        ]
        chrom_sizes = {"JAC": 7}
        result = RegionMapper(base_records, chrom_sizes).get_map()
        expected = [SeqRegion(base_pos=1, inf_pos=1, length=7)]
        self.assertEqual(expected, result["JAC"])

    def test_AltLongerThanRef_CorrectRegion(self):
        # base sequence:      T TAT    CGG
        # derived sequence:   T GCCAC  CGG
        base_records = [_MockVcfRecord(pos=2, ref="TAT", alts=["GCCAC"])]
        chrom_sizes = {"JAC": 7}

        result = RegionMapper(base_records, chrom_sizes).get_map()

        expected = [
            SeqRegion(base_pos=1, inf_pos=1, length=1),
            SeqRegion(
                base_pos=2,
                inf_pos=2,
                length=5,
                vcf_record_ref="TAT",
                vcf_record_alt="GCCAC",
            ),
            SeqRegion(base_pos=5, inf_pos=7, length=3),
        ]
        self.assertEqual(expected, result["JAC"])

    def test_TwoRecords_CorrectRegions(self):
        # base sequence:      T TAT    C G   G
        # derived sequence:   T GCCAC  C TTT G
        base_records = [
            _MockVcfRecord(pos=2, ref="TAT", alts=["GCCAC"]),
            _MockVcfRecord(pos=6, ref="G", alts=["TTT"]),
        ]

        chrom_sizes = {"JAC": 7}
        result = RegionMapper(base_records, chrom_sizes).get_map()

        expected = [
            SeqRegion(base_pos=1, inf_pos=1, length=1),
            SeqRegion(
                base_pos=2,
                inf_pos=2,
                length=5,
                vcf_record_ref="TAT",
                vcf_record_alt="GCCAC",
            ),
            SeqRegion(base_pos=5, inf_pos=7, length=1),
            SeqRegion(
                base_pos=6,
                inf_pos=8,
                length=3,
                vcf_record_ref="G",
                vcf_record_alt="TTT",
            ),
            SeqRegion(base_pos=7, inf_pos=11, length=1),
        ]

        self.assertEqual(expected, result["JAC"])

    def test_ThreeAdjacentRecords_CorrectRegions(self):
        # base sequence:      T TAT    C   G  G
        # derived sequence:   T GCCAC  TCT AA G
        base_records = [
            _MockVcfRecord(pos=2, ref="TAT", alts=["GCCAC"]),
            _MockVcfRecord(pos=5, ref="C", alts=["TCT"]),
            _MockVcfRecord(pos=6, ref="G", alts=["AA"]),
        ]
        chrom_sizes = {"JAC": 7}
        result = RegionMapper(base_records, chrom_sizes).get_map()

        expected = [
            SeqRegion(base_pos=1, inf_pos=1, length=1),
            SeqRegion(
                base_pos=2,
                inf_pos=2,
                length=5,
                vcf_record_ref="TAT",
                vcf_record_alt="GCCAC",
            ),
            SeqRegion(
                base_pos=5,
                inf_pos=7,
                length=3,
                vcf_record_ref="C",
                vcf_record_alt="TCT",
            ),
            SeqRegion(
                base_pos=6,
                inf_pos=10,
                length=2,
                vcf_record_ref="G",
                vcf_record_alt="AA",
            ),
            SeqRegion(base_pos=7, inf_pos=12, length=1),
        ]

        self.assertEqual(expected, list(result.values())[0])

    def test_TwoRecords_TwoDifferentChroms(self):
        # CHROM 1:
        #   base sequence:      GAA ATTC CAA
        #   derived sequence:   GAA A    CAA

        # CHROM 2:
        #   base sequence:      GCGCA A   CG
        #   derived sequence:   GCGCA AAC CG

        base_records = [
            _MockVcfRecord(pos=4, ref="ATTC", alts=["A"], chrom="Chrom_1"),
            _MockVcfRecord(pos=6, ref="A", alts=["AAC"], chrom="Chrom_2"),
        ]

        chrom_sizes = {"Chrom_1": 10, "Chrom_2": 8}
        result = RegionMapper(base_records, chrom_sizes).get_map()

        expected_Chrom_1 = [
            SeqRegion(base_pos=1, inf_pos=1, length=3),
            SeqRegion(
                base_pos=4,
                inf_pos=4,
                length=1,
                vcf_record_ref="ATTC",
                vcf_record_alt="A",
            ),
            SeqRegion(base_pos=8, inf_pos=5, length=3),
        ]

        expected_Chrom_2 = [
            SeqRegion(base_pos=1, inf_pos=1, length=5),
            SeqRegion(
                base_pos=6,
                inf_pos=6,
                length=3,
                vcf_record_ref="A",
                vcf_record_alt="AAC",
            ),
            SeqRegion(base_pos=7, inf_pos=9, length=2),
        ]
        expectations = {"Chrom_1": expected_Chrom_1, "Chrom_2": expected_Chrom_2}
        for key in expectations:
            self.assertEqual(expectations[key], result[key])

    def test_NoRecords(self):
        """
        Imagining a vcf from `infer` has no records, `build` would not have succeeded in
        the first place, having not built a prg given no variants.
        """
        # base sequence:      TTATCGG
        # derived sequence:   TTATCGG
        chrom_sizes = {}
        base_records = []
        with self.assertRaises(ValueError):
            result = RegionMapper(base_records, chrom_sizes).get_map()

    def test_chrom_with_no_records(self):
        """
        Can find variants against a chromosome with no initial variation
        """
        # CHROM1:
        #   base sequence:      CAAC
        #   derived sequence:   CAAC
        # CHROM2:
        #   base sequence:      A T CAA
        #   derived sequence:   A A CAA
        base_records = [_MockVcfRecord(pos=2, ref="T", alts=["A"], chrom="Chrom_2")]

        chrom_sizes = {"Chrom_1": 4, "Chrom_2": 5}
        result = RegionMapper(base_records, chrom_sizes).get_map()

        expected_Chrom_1 = [SeqRegion(base_pos=1, inf_pos=1, length=4)]
        expected_Chrom_2 = [
            SeqRegion(base_pos=1, inf_pos=1, length=1),
            SeqRegion(
                base_pos=2, inf_pos=2, length=1, vcf_record_ref="T", vcf_record_alt="A"
            ),
            SeqRegion(base_pos=3, inf_pos=3, length=3),
        ]

        expectations = {"Chrom_1": expected_Chrom_1, "Chrom_2": expected_Chrom_2}
        for key in expectations:
            self.assertEqual(expectations[key], result[key])


class TestSearchRegions(TestCase):
    def test_RecordStartsAtRegionMid(self):
        """
        The record is also the very last record in the regions.
        """
        # base sequence:      T TAT CGG
        # derived sequence:   T G   CGG

        mapped_regions = [
            SeqRegion(base_pos=1, inf_pos=1, length=1),
            SeqRegion(
                base_pos=2,
                inf_pos=2,
                length=1,
                vcf_record_ref="TAT",
                vcf_record_alt="G",
            ),
            SeqRegion(base_pos=5, inf_pos=3, length=3),
        ]

        vcf_record = _MockVcfRecord(pos=4, ref="G", alts=["A"])
        result = _find_start_region_index(vcf_record, mapped_regions)

        expected = 2
        self.assertEqual(expected, result)

    def test_RecordStartsAtRegionStart(self):
        # base sequence:      T TAT    CGG
        # derived sequence:   T GCCAC  CGG

        mapped_regions = [
            SeqRegion(base_pos=1, inf_pos=1, length=1),
            SeqRegion(
                base_pos=2,
                inf_pos=2,
                length=5,
                vcf_record_ref="TAT",
                vcf_record_alt="GCCAC",
            ),
            SeqRegion(base_pos=5, inf_pos=7, length=3),
        ]

        vcf_record = _MockVcfRecord(pos=2, ref="GC", alts=["GA"])
        result = _find_start_region_index(vcf_record, mapped_regions)

        expected = 1
        self.assertEqual(expected, result)

    def test_RecordStartsAtRegionEnd(self):
        # base sequence:      T TAT    C G   G
        # derived sequence:   T GCCAC  C TTT G

        mapped_regions = [
            SeqRegion(base_pos=1, inf_pos=1, length=1),
            SeqRegion(
                base_pos=2,
                inf_pos=2,
                length=5,
                vcf_record_ref="TAT",
                vcf_record_alt="GCCAC",
            ),
            SeqRegion(base_pos=5, inf_pos=7, length=1),
            SeqRegion(
                base_pos=6,
                inf_pos=8,
                length=3,
                vcf_record_ref="G",
                vcf_record_alt="TTT",
            ),
            SeqRegion(base_pos=7, inf_pos=11, length=1),
        ]

        vcf_record = _MockVcfRecord(pos=10, ref="T", alts=["A"])
        result = _find_start_region_index(vcf_record, mapped_regions)

        expected = 3
        self.assertEqual(expected, result)
