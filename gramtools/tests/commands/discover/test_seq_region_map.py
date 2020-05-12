from unittest import TestCase
from tempfile import mkdtemp
from pathlib import Path
from shutil import rmtree

from gramtools.tests.mocks import _MockVcfRecord
from gramtools.commands.discover.seq_region_map import (
    SeqRegion,
    SearchableSeqRegionsMap,
    BisectTarget,
)
from gramtools.commands.discover.discover import SeqRegionMapper


class TestRegionMapping(TestCase):
    """
    Note: as stated in the original code, the length of a region is the length of the alt string if we have a vcf record.
    """

    def test_SingleBaseAlt_CorrectRegion(self):
        # base sequence:      T TAT CGG
        # derived sequence:   T G   CGG
        base_records = [_MockVcfRecord(pos=2, ref="TAT", alts=["G"])]

        chrom_sizes = {"JAC": 7}

        result = SeqRegionMapper(base_records, chrom_sizes).get_map()

        expected = [
            SeqRegion(base_ref_start=1, pers_ref_start=1, length=1),
            SeqRegion(
                base_ref_start=2,
                pers_ref_start=2,
                length=1,
                vcf_record_ref="TAT",
                vcf_record_alt="G",
            ),
            SeqRegion(base_ref_start=5, pers_ref_start=3, length=3),
        ]
        self.assertEqual(expected, result["JAC"])

    def test_ref_call_produces_invariant_region_only(self):
        # base sequence:      T TAT CGG
        # derived sequence:   ^^^^^^^^^
        base_records = [
            _MockVcfRecord(pos=2, ref="TAT", alts=["G"], samples=[{"GT": [0]}])
        ]
        chrom_sizes = {"JAC": 7}
        result = SeqRegionMapper(base_records, chrom_sizes).get_map()
        expected = [SeqRegion(base_ref_start=1, pers_ref_start=1, length=7)]
        self.assertEqual(expected, result["JAC"])

    def test_AltLongerThanRef_CorrectRegion(self):
        # base sequence:      T TAT    CGG
        # derived sequence:   T GCCAC  CGG
        base_records = [_MockVcfRecord(pos=2, ref="TAT", alts=["GCCAC"])]
        chrom_sizes = {"JAC": 7}

        result = SeqRegionMapper(base_records, chrom_sizes).get_map()

        expected = [
            SeqRegion(base_ref_start=1, pers_ref_start=1, length=1),
            SeqRegion(
                base_ref_start=2,
                pers_ref_start=2,
                length=5,
                vcf_record_ref="TAT",
                vcf_record_alt="GCCAC",
            ),
            SeqRegion(base_ref_start=5, pers_ref_start=7, length=3),
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
        result = SeqRegionMapper(base_records, chrom_sizes).get_map()

        expected = [
            SeqRegion(base_ref_start=1, pers_ref_start=1, length=1),
            SeqRegion(
                base_ref_start=2,
                pers_ref_start=2,
                length=5,
                vcf_record_ref="TAT",
                vcf_record_alt="GCCAC",
            ),
            SeqRegion(base_ref_start=5, pers_ref_start=7, length=1),
            SeqRegion(
                base_ref_start=6,
                pers_ref_start=8,
                length=3,
                vcf_record_ref="G",
                vcf_record_alt="TTT",
            ),
            SeqRegion(base_ref_start=7, pers_ref_start=11, length=1),
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
        result = SeqRegionMapper(base_records, chrom_sizes).get_map()

        expected = [
            SeqRegion(base_ref_start=1, pers_ref_start=1, length=1),
            SeqRegion(
                base_ref_start=2,
                pers_ref_start=2,
                length=5,
                vcf_record_ref="TAT",
                vcf_record_alt="GCCAC",
            ),
            SeqRegion(
                base_ref_start=5,
                pers_ref_start=7,
                length=3,
                vcf_record_ref="C",
                vcf_record_alt="TCT",
            ),
            SeqRegion(
                base_ref_start=6,
                pers_ref_start=10,
                length=2,
                vcf_record_ref="G",
                vcf_record_alt="AA",
            ),
            SeqRegion(base_ref_start=7, pers_ref_start=12, length=1),
        ]

        self.assertEqual(expected, list(result.values())[0])

    def test_TwoRecords_TwoDifferentChroms(self):
        base_records = [
            _MockVcfRecord(pos=4, ref="ATTC", alts=["A"], chrom="Chrom_1"),
            _MockVcfRecord(pos=6, ref="A", alts=["AAC"], chrom="Chrom_2"),
        ]

        chrom_sizes = {"Chrom_1": 10, "Chrom_2": 8}
        result = SeqRegionMapper(base_records, chrom_sizes).get_map()

        expected_Chrom_1 = [
            SeqRegion(base_ref_start=1, pers_ref_start=1, length=3),
            SeqRegion(
                base_ref_start=4,
                pers_ref_start=4,
                length=1,
                vcf_record_ref="ATTC",
                vcf_record_alt="A",
            ),
            SeqRegion(base_ref_start=8, pers_ref_start=5, length=3),
        ]

        expected_Chrom_2 = [
            SeqRegion(base_ref_start=1, pers_ref_start=1, length=5),
            SeqRegion(
                base_ref_start=6,
                pers_ref_start=6,
                length=3,
                vcf_record_ref="A",
                vcf_record_alt="AAC",
            ),
            SeqRegion(base_ref_start=7, pers_ref_start=9, length=2),
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
            result = SeqRegionMapper(base_records, chrom_sizes).get_map()

    def test_chrom_with_no_records(self):
        """
        Need to map chroms with no initial variation too
        """
        base_records = [_MockVcfRecord(pos=2, ref="T", alts=["A"], chrom="Chrom_2")]

        chrom_sizes = {"Chrom_1": 4, "Chrom_2": 5}
        result = SeqRegionMapper(base_records, chrom_sizes).get_map()

        expected_Chrom_1 = [SeqRegion(base_ref_start=1, pers_ref_start=1, length=4)]
        expected_Chrom_2 = [
            SeqRegion(base_ref_start=1, pers_ref_start=1, length=1),
            SeqRegion(
                base_ref_start=2,
                pers_ref_start=2,
                length=1,
                vcf_record_ref="T",
                vcf_record_alt="A",
            ),
            SeqRegion(base_ref_start=3, pers_ref_start=3, length=3),
        ]

        expectations = {"Chrom_1": expected_Chrom_1, "Chrom_2": expected_Chrom_2}
        for key in expectations:
            self.assertEqual(expectations[key], result[key])


class TestSearchRegions(TestCase):
    def test_base_ref_pers_ref_same_results(self):
        mapped_regions = [
            SeqRegion(base_ref_start=1, pers_ref_start=1, length=1),
            SeqRegion(
                base_ref_start=2,
                pers_ref_start=2,
                length=3,
                vcf_record_ref="TAT",
                vcf_record_alt="GCC",
            ),
            SeqRegion(base_ref_start=5, pers_ref_start=5, length=3),
        ]

        vcf_record_in_var_region = _MockVcfRecord(pos=2, ref="GC", alts=["GA"])
        vcf_record_in_nonvar_region = _MockVcfRecord(pos=1, ref="A", alts=["T"])
        searcher = SearchableSeqRegionsMap({"JAC": mapped_regions})

        for target in BisectTarget:
            result = searcher.bisect("JAC", vcf_record_in_var_region.pos, target)
            self.assertEqual(1, result)

        for target in BisectTarget:
            result = searcher.bisect("JAC", vcf_record_in_nonvar_region.pos, target)
            self.assertEqual(0, result)

    def test_pers_ref_further_than_base_ref(self):
        mapped_regions = [
            SeqRegion(base_ref_start=1, pers_ref_start=1, length=1),
            SeqRegion(
                base_ref_start=2,
                pers_ref_start=2,
                length=1,
                vcf_record_ref="TAT",
                vcf_record_alt="G",
            ),
            SeqRegion(base_ref_start=5, pers_ref_start=3, length=3),
        ]

        vcf_record = _MockVcfRecord(pos=4, ref="G", alts=["A"])
        searcher = SearchableSeqRegionsMap({"JAC": mapped_regions})

        # Pers ref bisection
        pers_ref_result = searcher.bisect("JAC", vcf_record.pos, BisectTarget.PERS_REF)
        self.assertEqual(2, pers_ref_result)

        # Base ref bisection
        base_ref_result = searcher.bisect("JAC", vcf_record.pos, BisectTarget.BASE_REF)
        self.assertEqual(1, base_ref_result)

    def test_base_ref_further_than_pers_ref(self):
        mapped_regions = [
            SeqRegion(base_ref_start=1, pers_ref_start=1, length=1),
            SeqRegion(
                base_ref_start=2,
                pers_ref_start=2,
                length=5,
                vcf_record_ref="TAT",
                vcf_record_alt="GCCAC",
            ),
            SeqRegion(
                base_ref_start=5,
                pers_ref_start=7,
                length=3,
                vcf_record_ref="G",
                vcf_record_alt="TTT",
            ),
        ]

        vcf_record = _MockVcfRecord(pos=6, ref="T", alts=["A"])
        searcher = SearchableSeqRegionsMap({"JAC": mapped_regions})

        pers_ref_result = searcher.bisect("JAC", vcf_record.pos, BisectTarget.PERS_REF)
        self.assertEqual(1, pers_ref_result)

        base_ref_result = searcher.bisect("JAC", vcf_record.pos, BisectTarget.BASE_REF)
        self.assertEqual(2, base_ref_result)

    def test_retrieve_searched_region(self):
        mapped_regions_1 = [
            SeqRegion(base_ref_start=1, pers_ref_start=1, length=1),
            SeqRegion(
                base_ref_start=2,
                pers_ref_start=2,
                length=5,
                vcf_record_ref="TAT",
                vcf_record_alt="GCCAC",
            ),
        ]

        mapped_regions_2 = [SeqRegion(base_ref_start=1, pers_ref_start=1, length=200)]

        searcher = SearchableSeqRegionsMap(
            {"chr1": mapped_regions_1, "chr2": mapped_regions_2}
        )

        vcf_record = _MockVcfRecord(pos=100, ref="T", alts=["A"], chrom="chr2")
        its_index = searcher.bisect(
            vcf_record.chrom, vcf_record.pos, BisectTarget.PERS_REF
        )
        self.assertEqual(
            searcher.get_region(vcf_record.chrom, its_index), mapped_regions_2[0]
        )

        vcf_record = _MockVcfRecord(pos=4, ref="C", alts=["A"], chrom="chr1")
        its_index = searcher.bisect(
            vcf_record.chrom, vcf_record.pos, BisectTarget.PERS_REF
        )
        self.assertEqual(
            searcher.get_region(vcf_record.chrom, its_index), mapped_regions_1[1]
        )


class TestSearcherSerialisation(TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.mapped_regions = [
            SeqRegion(base_ref_start=1, pers_ref_start=1, length=1),
            SeqRegion(
                base_ref_start=2,
                pers_ref_start=2,
                length=2,
                vcf_record_ref="TAT",
                vcf_record_alt="CC",
            ),
        ]

    def test_json_serialise_SeqRegions(self):
        reg1 = self.mapped_regions[0].to_json()
        expected = {
            "SeqRegion": {"base_ref_start": 1, "pers_ref_start": 1, "length": 1}
        }
        self.assertEqual(expected, reg1)

        reg2 = self.mapped_regions[1].to_json()
        expected = {
            "SeqRegion": {
                "base_ref_start": 2,
                "pers_ref_start": 2,
                "length": 2,
                "vcf_record_ref": "TAT",
                "vcf_record_alt": "CC",
            }
        }
        self.assertEqual(expected, reg2)

    def test_json_deserialise_SeqRegion(self):
        reg1 = {"SeqRegion": {"base_ref_start": 1, "pers_ref_start": 1, "length": 1}}
        self.assertEqual(self.mapped_regions[0], SeqRegion.from_json(reg1))

        reg2 = {
            "SeqRegion": {
                "base_ref_start": 2,
                "pers_ref_start": 2,
                "length": 2,
                "vcf_record_ref": "TAT",
                "vcf_record_alt": "CC",
            }
        }
        self.assertEqual(self.mapped_regions[1], SeqRegion.from_json(reg2))

    def test_dump_and_load_recapitulates_map(self):
        searcher = SearchableSeqRegionsMap({"JAC": self.mapped_regions})
        tmpdir = Path(mkdtemp())
        tmpfile = tmpdir / "map.json"
        searcher.dump_to(tmpfile)
        loaded_searcher = SearchableSeqRegionsMap.load_from(tmpfile)
        self.assertEqual(searcher, loaded_searcher)
        rmtree(tmpdir)

    def test_dump_and_load_without_sequences(self):
        """
        Serialisation without the REF and ALT SeqRegion sequences
        """
        searcher = SearchableSeqRegionsMap({"JAC": self.mapped_regions})
        tmpdir = Path(mkdtemp())
        tmpfile = tmpdir / "map.json"
        searcher.dump_to(tmpfile, dump_sequences=False)
        loaded_searcher = SearchableSeqRegionsMap.load_from(tmpfile)

        self.assertEqual(
            searcher.get_region("JAC", 0), loaded_searcher.get_region("JAC", 0)
        )
        self.assertEqual(SeqRegion(2, 2, 2), loaded_searcher.get_region("JAC", 1))
        rmtree(tmpdir)
