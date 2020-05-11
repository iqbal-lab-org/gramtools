from unittest import TestCase, main as unittest_main

from gramtools.tests.mocks import _MockVcfRecord
from gramtools.commands.discover import discover
from gramtools.commands.discover.seq_region_map import SearchableSeqRegionsMap


class TestRebaseVcfRecord(TestCase):
    def test_SingleSNPInNonSite(self):
        # base sequence:      T TAT CGG
        # derived sequence:   T G   CGG
        chrom_sizes = {"JAC": 5}
        base_records = [_MockVcfRecord(pos=2, ref="TAT", alts=["G"])]
        region_map = discover.SeqRegionMapper(base_records, chrom_sizes).get_map()
        region_searcher = SearchableSeqRegionsMap(region_map)

        discov_record = _MockVcfRecord(pos=3, ref="C", alts=["G"])
        new_vcf_record = discover._rebase_vcf_record(
            discov_record, "JAC", region_searcher
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
        mapped_regions = discover.SeqRegionMapper(base_records, chrom_sizes).get_map()
        region_searcher = SearchableSeqRegionsMap(mapped_regions)

        discov_record = _MockVcfRecord(pos=1, ref="TT", alts=["GA"], chrom="chr2")
        new_vcf_record = discover._rebase_vcf_record(
            discov_record, "chr2", region_searcher
        )
        self.assertEqual(discov_record, new_vcf_record)

    def test_StartsAtNonSite_EndsAtSite(self):
        # base sequence:      T TAT CGG
        # derived sequence:   T G   CGG
        chrom_sizes = {"JAC": 7}
        base_records = [_MockVcfRecord(pos=2, ref="TAT", alts=["G"])]
        mapped_regions = discover.SeqRegionMapper(base_records, chrom_sizes).get_map()
        region_searcher = SearchableSeqRegionsMap(mapped_regions)

        discov_record = _MockVcfRecord(pos=1, ref="TG", alts=["TAA"])
        new_vcf_record = discover._rebase_vcf_record(
            discov_record, "JAC", region_searcher
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
        mapped_regions = discover.SeqRegionMapper(base_records, chrom_sizes).get_map()
        region_searcher = SearchableSeqRegionsMap(mapped_regions)

        discov_record = _MockVcfRecord(pos=1, ref="TGCG", alts=["GGCT"])

        new_vcf_record = discover._rebase_vcf_record(
            discov_record, "JAC", region_searcher
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
        mapped_regions = discover.SeqRegionMapper(base_records, chrom_sizes).get_map()
        region_searcher = SearchableSeqRegionsMap(mapped_regions)

        discov_record = _MockVcfRecord(pos=9, ref="G", alts=["A"])

        new_vcf_record = discover._rebase_vcf_record(
            discov_record, "JAC", region_searcher
        )

        result = _MockVcfRecord(
            new_vcf_record.pos, new_vcf_record.ref, new_vcf_record.alts
        )
        expected = _MockVcfRecord(8, "T", ["TCTAC"])

        self.assertEqual(expected, result)

    def test_multiple_deletions(self):
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
        mapped_regions = discover.SeqRegionMapper(base_records, chrom_sizes).get_map()
        region_searcher = SearchableSeqRegionsMap(mapped_regions)

        discov_record = _MockVcfRecord(pos=4, ref="ATC", alts=["A"])
        new_vcf_record = discover._rebase_vcf_record(
            discov_record, "JAC", region_searcher
        )

        result = _MockVcfRecord(
            new_vcf_record.pos, new_vcf_record.ref, new_vcf_record.alts
        )
        expected = _MockVcfRecord(pos=5, ref="GCTAC", alts=["GA"])

        self.assertEqual(expected, result)


if __name__ == "__main__":
    unittest_main()
