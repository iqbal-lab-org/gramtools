import unittest
import os
import collections

from gramtools.commands.discover import discover
from gramtools.commands.discover.region_mapper import _Region
from gramtools.utils import prg_local_parser
from gramtools.tests.mocks import _MockVcfRecord


class TestRegionMapping(unittest.TestCase):
    """
    Note: as stated in the original code, the length of a region is the length of the (first) alt string if we have a vcf record.
    """

    def test_SingleBaseAlt_CorrectRegion(self):
        # base sequence:      T TAT CGG
        # derived sequence:   T G   CGG
        base_records = [_MockVcfRecord(pos=2, ref="TAT", alts=["G"])]

        chrom_sizes = {"JAC": 7}

        result = discover.RegionMapper(base_records, chrom_sizes).get_mapped()

        expected = [
            _Region(base_pos=1, inf_pos=1, length=1),
            _Region(
                base_pos=2,
                inf_pos=2,
                length=1,
                vcf_record_ref="TAT",
                vcf_record_alt="G",
            ),
            _Region(base_pos=5, inf_pos=3, length=3),
        ]
        self.assertEqual(expected, list(result.values())[0])

    def test_AltLongerThanRef_CorrectRegion(self):
        # base sequence:      T TAT    CGG
        # derived sequence:   T GCCAC  CGG
        base_records = [_MockVcfRecord(pos=2, ref="TAT", alts=["GCCAC"])]
        chrom_sizes = {"JAC": 7}

        result = discover.RegionMapper(base_records, chrom_sizes).get_mapped()

        expected = [
            _Region(base_pos=1, inf_pos=1, length=1),
            _Region(
                base_pos=2,
                inf_pos=2,
                length=5,
                vcf_record_ref="TAT",
                vcf_record_alt="GCCAC",
            ),
            _Region(base_pos=5, inf_pos=7, length=3),
        ]
        self.assertEqual(expected, list(result.values())[0])

    def test_TwoRecords_CorrectRegions(self):
        # base sequence:      T TAT    C G   G
        # derived sequence:   T GCCAC  C TTT G
        base_records = [
            _MockVcfRecord(pos=2, ref="TAT", alts=["GCCAC"]),
            _MockVcfRecord(pos=6, ref="G", alts=["TTT"]),
        ]

        chrom_sizes = {"JAC": 7}
        result = discover.RegionMapper(base_records, chrom_sizes).get_mapped()

        expected = [
            _Region(base_pos=1, inf_pos=1, length=1),
            _Region(
                base_pos=2,
                inf_pos=2,
                length=5,
                vcf_record_ref="TAT",
                vcf_record_alt="GCCAC",
            ),
            _Region(base_pos=5, inf_pos=7, length=1),
            _Region(
                base_pos=6,
                inf_pos=8,
                length=3,
                vcf_record_ref="G",
                vcf_record_alt="TTT",
            ),
            _Region(base_pos=7, inf_pos=11, length=1),
        ]

        self.assertEqual(expected, list(result.values())[0])

    def test_ThreeAdjacentRecords_CorrectRegions(self):
        # base sequence:      T TAT    C   G  G
        # derived sequence:   T GCCAC  TCT AA G
        base_records = [
            _MockVcfRecord(pos=2, ref="TAT", alts=["GCCAC"]),
            _MockVcfRecord(pos=5, ref="C", alts=["TCT"]),
            _MockVcfRecord(pos=6, ref="G", alts=["AA"]),
        ]
        chrom_sizes = {"JAC": 7}
        result = discover.RegionMapper(base_records, chrom_sizes).get_mapped()

        expected = [
            _Region(base_pos=1, inf_pos=1, length=1),
            _Region(
                base_pos=2,
                inf_pos=2,
                length=5,
                vcf_record_ref="TAT",
                vcf_record_alt="GCCAC",
            ),
            _Region(
                base_pos=5,
                inf_pos=7,
                length=3,
                vcf_record_ref="C",
                vcf_record_alt="TCT",
            ),
            _Region(
                base_pos=6,
                inf_pos=10,
                length=2,
                vcf_record_ref="G",
                vcf_record_alt="AA",
            ),
            _Region(base_pos=7, inf_pos=12, length=1),
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
        result = discover.RegionMapper(base_records, chrom_sizes).get_mapped()

        expected_Chrom_1 = [
            _Region(base_pos=1, inf_pos=1, length=3),
            _Region(
                base_pos=4,
                inf_pos=4,
                length=1,
                vcf_record_ref="ATTC",
                vcf_record_alt="A",
            ),
            _Region(base_pos=8, inf_pos=5, length=3),
        ]

        expected_Chrom_2 = [
            _Region(base_pos=1, inf_pos=1, length=5),
            _Region(
                base_pos=6,
                inf_pos=6,
                length=3,
                vcf_record_ref="A",
                vcf_record_alt="AAC",
            ),
            _Region(base_pos=7, inf_pos=9, length=2),
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
            result = discover.RegionMapper(base_records, chrom_sizes).get_mapped()

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
        result = discover.RegionMapper(base_records, chrom_sizes).get_mapped()

        expected_Chrom_1 = [_Region(base_pos=1, inf_pos=1, length=4)]
        expected_Chrom_2 = [
            _Region(base_pos=1, inf_pos=1, length=1),
            _Region(
                base_pos=2, inf_pos=2, length=1, vcf_record_ref="T", vcf_record_alt="A"
            ),
            _Region(base_pos=3, inf_pos=3, length=3),
        ]

        expectations = {"Chrom_1": expected_Chrom_1, "Chrom_2": expected_Chrom_2}
        for key in expectations:
            self.assertEqual(expectations[key], result[key])


class TestSearchRegions(unittest.TestCase):
    def test_RecordStartsAtRegionMid(self):
        """
        The record is also the very last record in the regions.
        """
        # base sequence:      T TAT CGG
        # derived sequence:   T G   CGG
        base_records = [_MockVcfRecord(pos=2, ref="TAT", alts=["G"])]

        chrom_sizes = {"JAC": 7}

        mapped_regions = [
            _Region(base_pos=1, inf_pos=1, length=1),
            _Region(
                base_pos=2,
                inf_pos=2,
                length=1,
                vcf_record_ref="TAT",
                vcf_record_alt="G",
            ),
            _Region(base_pos=5, inf_pos=3, length=3),
        ]

        vcf_record = _MockVcfRecord(pos=4, ref="G", alts=["A"])
        result = discover._find_start_region_index(vcf_record, mapped_regions)

        expected = 2
        self.assertEqual(expected, result)

    def test_RecordStartsAtRegionStart(self):
        # base sequence:      T TAT    CGG
        # derived sequence:   T GCCAC  CGG
        base_records = [_MockVcfRecord(pos=2, ref="TAT", alts=["GCCAC"])]

        mapped_regions = [
            _Region(base_pos=1, inf_pos=1, length=1),
            _Region(
                base_pos=2,
                inf_pos=2,
                length=5,
                vcf_record_ref="TAT",
                vcf_record_alt="GCCAC",
            ),
            _Region(base_pos=5, inf_pos=7, length=3),
        ]

        vcf_record = _MockVcfRecord(pos=2, ref="GC", alts=["GA"])
        result = discover._find_start_region_index(vcf_record, mapped_regions)

        expected = 1
        self.assertEqual(expected, result)

    def test_RecordStartsAtRegionEnd(self):
        # base sequence:      T TAT    C G   G
        # derived sequence:   T GCCAC  C TTT G
        base_records = [
            _MockVcfRecord(pos=2, ref="TAT", alts=["GCCAC"]),
            _MockVcfRecord(pos=6, ref="G", alts=["TTT"]),
        ]

        mapped_regions = [
            _Region(base_pos=1, inf_pos=1, length=1),
            _Region(
                base_pos=2,
                inf_pos=2,
                length=5,
                vcf_record_ref="TAT",
                vcf_record_alt="GCCAC",
            ),
            _Region(base_pos=5, inf_pos=7, length=1),
            _Region(
                base_pos=6,
                inf_pos=8,
                length=3,
                vcf_record_ref="G",
                vcf_record_alt="TTT",
            ),
            _Region(base_pos=7, inf_pos=11, length=1),
        ]

        vcf_record = _MockVcfRecord(pos=10, ref="T", alts=["A"])
        result = discover._find_start_region_index(vcf_record, mapped_regions)

        expected = 3
        self.assertEqual(expected, result)


class TestRebaseVcfRecord(unittest.TestCase):
    def test_SingleSNPInNonSite(self):
        # base sequence:      T TAT CGG
        # derived sequence:   T G   CGG
        chrom_sizes = {"JAC": 5}
        base_records = [_MockVcfRecord(pos=2, ref="TAT", alts=["G"])]
        mapped_regions = discover.RegionMapper(base_records, chrom_sizes).get_mapped()

        discov_record = _MockVcfRecord(pos=3, ref="C", alts=["G"])
        new_vcf_record = discover._rebase_vcf_record(
            discov_record, list(mapped_regions.values())[0]
        )

        result = _MockVcfRecord(
            new_vcf_record.pos, new_vcf_record.ref, new_vcf_record.alts
        )
        expected = _MockVcfRecord(pos=5, ref="C", alts=["G"])

        self.assertEqual(expected, result)

    def test_variant_in_chromo_with_no_prg_variants(self):
        # chr1 base:    T TAT CGG
        # chr1 derived: T G   CGG
        # chr2 base:    TTTTT
        # chr2 derived: TTTTT

        chrom_sizes = {"chr1": 7, "chr2": 5}
        base_records = [_MockVcfRecord(pos=2, ref="TAT", alts=["G"], chrom="chr1")]
        mapped_regions = discover.RegionMapper(base_records, chrom_sizes).get_mapped()

        discov_record = _MockVcfRecord(pos=1, ref="TT", alts=["GA"], chrom="chr2")
        new_vcf_record = discover._rebase_vcf_record(
            discov_record, mapped_regions["chr2"]
        )
        self.assertEqual(discov_record, discov_record)

    def test_StartsAtNonSite_EndsAtSite(self):
        # base sequence:      T TAT CGG
        # derived sequence:   T G   CGG
        chrom_sizes = {"JAC": 7}
        base_records = [_MockVcfRecord(pos=2, ref="TAT", alts=["G"])]
        mapped_regions = discover.RegionMapper(base_records, chrom_sizes).get_mapped()

        discov_record = _MockVcfRecord(pos=1, ref="TG", alts=["TAA"])
        new_vcf_record = discover._rebase_vcf_record(
            discov_record, mapped_regions["JAC"]
        )

        result = _MockVcfRecord(
            new_vcf_record.pos, new_vcf_record.ref, new_vcf_record.alts
        )
        expected = _MockVcfRecord(1, "TTAT", ["TAA"])

        self.assertEqual(expected, result)

    def test_SiteInBetweenNonSites(self):
        """
        A test case where the variation on top of the inferred reference overlaps: a non-variant site, a variant site,
        and a non-variant site in the prg.

        What we need is for the rebased ref to include all three sites.
        """
        # base sequ: T TAT CGG
        # secondary: T G   CGG
        chrom_sizes = {"JAC": 7}
        base_records = [_MockVcfRecord(pos=2, ref="TAT", alts=["G"])]

        discov_record = _MockVcfRecord(pos=1, ref="TGCG", alts=["GGCT"])

        mapped_regions = discover.RegionMapper(base_records, chrom_sizes).get_mapped()
        new_vcf_record = discover._rebase_vcf_record(
            discov_record, list(mapped_regions.values())[0]
        )

        result = _MockVcfRecord(
            new_vcf_record.pos, new_vcf_record.ref, new_vcf_record.alts
        )
        expected = _MockVcfRecord(pos=1, ref="TTATCG", alts=["GGCT"])

        self.assertEqual(expected, result)

    def test_SNP_OnTopOfIndel(self):
        """
        A test case where we find a SNP on top of an insertion in the inferred reference.

        What we need is for the rebased alt to include the flanking alt bases, which are implied to be present in the discov_record.
        """
        # base sequ: T TAT CGG T     A
        # secondary: T G   CGG TCTGC A
        chrom_sizes = {"JAC": 9}
        base_records = [
            _MockVcfRecord(pos=2, ref="TAT", alts=["G"]),
            _MockVcfRecord(pos=8, ref="T", alts=["TCTGC"]),
        ]

        discov_record = _MockVcfRecord(pos=9, ref="G", alts=["A"])

        mapped_regions = discover.RegionMapper(base_records, chrom_sizes).get_mapped()
        new_vcf_record = discover._rebase_vcf_record(
            discov_record, list(mapped_regions.values())[0]
        )

        result = _MockVcfRecord(
            new_vcf_record.pos, new_vcf_record.ref, new_vcf_record.alts
        )
        expected = _MockVcfRecord(8, "T", ["TCTAC"])

        self.assertEqual(expected, result)

    def test_Deletion_OnTopOfDeletion_WithExtraDeletionInNonSite(self):
        """
        A test case where we discover a deletion on top of a deletion in a variant site;
        as well as an extra deletion in a non-variant site.

        There is also a SNP among the original deletion, to make it plausible that quasimap/infer picks this variant.

        To make it harder, the discovered variation is also reported inside a variant site, so we expect the rebased alt to be elongated.

        We expect the rebased ref to include all deleted bases.
        """
        # base reference:     CAA C GCTA CAA
        # inferred reference: C   C GAT  CAA

        chrom_sizes = {"JAC": 11}
        base_records = [
            _MockVcfRecord(pos=1, ref="CAA", alts=["C"]),
            _MockVcfRecord(pos=5, ref="GCTA", alts=["GAT"]),
        ]

        discov_record = _MockVcfRecord(pos=4, ref="ATC", alts=["A"])

        mapped_regions = discover.RegionMapper(base_records, chrom_sizes).get_mapped()
        new_vcf_record = discover._rebase_vcf_record(
            discov_record, list(mapped_regions.values())[0]
        )

        result = _MockVcfRecord(
            new_vcf_record.pos, new_vcf_record.ref, new_vcf_record.alts
        )
        expected = _MockVcfRecord(pos=5, ref="GCTAC", alts=["GA"])

        self.assertEqual(expected, result)


class TestGetBaseReference(unittest.TestCase):
    """
    Pick the first allele of each variant site of the prg to make 'base' reference.
    """

    def test_GivenPrgWithTwoSites_CorrectReferenceExtracted(self):
        test_prg = "test_prg"
        test_inferred = "test_inferred"

        prg_seq = "TT5A6T5AA7C8A7"

        allele_indexes = (0 for _ in range(4))

        with open(test_prg, "w") as prg:
            prg.write(prg_seq)

        prg_parser = prg_local_parser.Prg_Local_Parser(
            test_prg, test_inferred, fasta_header="", allele_indexes=allele_indexes
        )
        prg_parser.parse()

        with open(test_inferred, "r") as inferred:
            result = inferred.readlines()[-1].strip()

        expected = "TTAAAC"
        os.remove(test_prg)
        os.remove(test_inferred)

        self.assertEqual(expected, result)


if __name__ == "__main__":
    unittest.main()
