import unittest
import vcf
import os

from ...commands import discover
from ... import prg_local_parser


class TestSecondaryRegions(unittest.TestCase):
    """
    Note: as stated in the original code, the length of a region is the length of the (first) ALT string if we have a vcf record.
    """
    def test_SingleBaseAlt_CorrectRegion(self):
        # base sequence:      T TAT CGG
        # secondary sequence: T G   CGG
        base_records = [
            _MockVcfRecord(POS=2, REF="TAT", ALT=['G'])
        ]

        secondary_sequence_length = 5

        result = discover._flag_personalised_reference_regions(base_records, secondary_sequence_length)

        expected = [
            discover._Region(base_POS = 1, inf_POS = 1, length = 1),
            discover._Region(base_POS = 2, inf_POS = 2, length = 1, vcf_record_REF = 'TAT', vcf_record_ALT = 'G'),
            discover._Region(base_POS = 5, inf_POS = 3, length = 3)
        ]
        self.assertEqual(expected, result)


    def test_AltLongerThanRef_CorrectRegion(self):
        # base sequence:      T TAT    CGG
        # secondary sequence: T GCCAC  CGG
        base_records = [
            _MockVcfRecord(POS=2, REF="TAT", ALT=['GCCAC'])
        ]
        secondary_sequence_length = 9

        result = discover._flag_personalised_reference_regions(base_records, secondary_sequence_length)

        expected = [
            discover._Region(base_POS = 1, inf_POS = 1, length = 1),
            discover._Region(base_POS = 2, inf_POS = 2, length = 5, vcf_record_REF = 'TAT', vcf_record_ALT = 'GCCAC'),
            discover._Region(base_POS = 5, inf_POS = 7, length = 3)
        ]
        self.assertEqual(expected, result)


    def test_TwoRecords_CorrectRegions(self):
        # base sequence:      T TAT    C G   G
        # secondary sequence: T GCCAC  C TTT G
        base_records = [
            _MockVcfRecord(POS = 2, REF  ="TAT", ALT = ["GCCAC"]),
            _MockVcfRecord(POS = 6, REF = "G", ALT = ["TTT"])
        ]

        secondary_reference_length = 11
        result = discover._flag_personalised_reference_regions(base_records, secondary_reference_length)

        expected = [
            discover._Region(base_POS = 1, inf_POS = 1, length = 1),
            discover._Region(base_POS = 2, inf_POS = 2, length = 5, vcf_record_REF = "TAT", vcf_record_ALT = 'GCCAC'),
            discover._Region(base_POS = 5, inf_POS = 7, length = 1),
            discover._Region(base_POS = 6, inf_POS = 8, length = 3, vcf_record_REF = "G", vcf_record_ALT = 'TTT'),
            discover._Region(base_POS = 7, inf_POS = 11, length = 1)
        ]

        self.assertEqual(expected, result)


    def test_ThreeAdjacentRecords_CorrectRegions(self):
        # base sequence:      T TAT    C   G  G
        # secondary sequence: T GCCAC  TCT AA G
        base_records = [
            _MockVcfRecord(POS=2, REF="TAT", ALT=['GCCAC']),
            _MockVcfRecord(POS=5, REF="C", ALT=['TCT']),
            _MockVcfRecord(POS=6, REF="G", ALT=['AA']),
        ]
        secondary_reference_length = 12
        result = discover._flag_personalised_reference_regions(base_records, secondary_reference_length)

        expected = [
            discover._Region(base_POS = 1, inf_POS = 1, length = 1),
            discover._Region(base_POS = 2, inf_POS = 2, length = 5, vcf_record_REF = "TAT", vcf_record_ALT = 'GCCAC'),
            discover._Region(base_POS = 5, inf_POS = 7, length = 3, vcf_record_REF = "C", vcf_record_ALT = "TCT"),
            discover._Region(base_POS = 6, inf_POS = 10, length = 2, vcf_record_REF = "G", vcf_record_ALT = 'AA'),
            discover._Region(base_POS = 7, inf_POS = 12, length = 1)
        ]

        self.assertEqual(expected, result)


    def test_NoRecords(self):
        # base sequence:      TTATCGG
        # secondary sequence: TTATCGG
        secondary_reference_length = 7
        base_records = []
        result = discover._flag_personalised_reference_regions(base_records, secondary_reference_length)

        expected = [discover._Region(base_POS = 1, inf_POS = 1, length = 7)]
        self.assertEqual(expected, result)


class TestSearchRegions(unittest.TestCase):
    """
    Note: the secondary_regions spelled out are tested for correctness in upstream tests: see @class TestSecondaryRegions
    So, we can safely use them if those tests don't break.

    Note(2): we find the first matching region by matching the vcf_record's POS with the POS coordinates of the inferred
    reference.
    """
    def test_RecordStartsAtRegionMid(self):
        """
        The record is also the very last record in the regions.
        """
        # base sequence:      T TAT CGG
        # secondary sequence: T G   CGG
        base_records = [
            _MockVcfRecord(POS=2, REF="TAT", ALT=['G'])
        ]

        secondary_sequence_length = 5

        secondary_regions = [
            discover._Region(base_POS = 1, inf_POS = 1, length = 1),
            discover._Region(base_POS = 2, inf_POS = 2, length = 1, vcf_record_REF = 'TAT', vcf_record_ALT = 'G'),
            discover._Region(base_POS = 5, inf_POS = 3, length = 3)
        ]

        vcf_record = _MockVcfRecord(POS=4, REF='G', ALT=['A'])
        result = discover._find_start_region_index(vcf_record, secondary_regions)

        expected = 2
        self.assertEqual(expected, result)

    def test_RecordStartsAtRegionStart(self):
        # base sequence:      T TAT    CGG
        # secondary sequence: T GCCAC  CGG
        base_records = [
            _MockVcfRecord(POS=2, REF="TAT", ALT=['GCCAC'])
        ]

        secondary_regions = [
            discover._Region(base_POS = 1, inf_POS = 1, length = 1),
            discover._Region(base_POS = 2, inf_POS = 2, length = 5, vcf_record_REF = 'TAT', vcf_record_ALT = 'GCCAC'),
            discover._Region(base_POS = 5, inf_POS = 7, length = 3)
        ]

        vcf_record = _MockVcfRecord(POS=2, REF='GC', ALT = ['GA'])
        result = discover._find_start_region_index(vcf_record, secondary_regions)

        expected = 1
        self.assertEqual(expected, result)


    def test_RecordStartsAtRegionEnd(self):
        # base sequence:      T TAT    C G   G
        # secondary sequence: T GCCAC  C TTT G
        base_records = [
            _MockVcfRecord(POS = 2, REF  ="TAT", ALT = ["GCCAC"]),
            _MockVcfRecord(POS = 6, REF = "G", ALT = ["TTT"])
        ]

        secondary_regions = [
            discover._Region(base_POS = 1, inf_POS = 1, length = 1),
            discover._Region(base_POS = 2, inf_POS = 2, length = 5, vcf_record_REF = "TAT", vcf_record_ALT = 'GCCAC'),
            discover._Region(base_POS = 5, inf_POS = 7, length = 1),
            discover._Region(base_POS = 6, inf_POS = 8, length = 3, vcf_record_REF = "G", vcf_record_ALT = 'TTT'),
            discover._Region(base_POS = 7, inf_POS = 11, length = 1)
        ]

        vcf_record = _MockVcfRecord(POS=10, REF='T', ALT=['A'])
        result = discover._find_start_region_index(vcf_record, secondary_regions)

        expected = 3
        self.assertEqual(expected, result)




class TestRebaseVcfRecord(unittest.TestCase):
    def test_SingleSNPInNonSite(self):
        # base sequence:      T TAT CGG
        # secondary sequence: T G   CGG
        secondary_reference_length = 5
        base_records = [
            _MockVcfRecord(POS=2, REF="TAT", ALT=['G'])
        ]
        secondary_regions = discover._flag_personalised_reference_regions(base_records, secondary_reference_length)

        secondary_vcf_record = _MockVcfRecord(POS=3, REF='C', ALT=['G'])
        new_vcf_record = discover._rebase_vcf_record(secondary_vcf_record, secondary_regions)

        result = _MockVcfRecord(new_vcf_record.POS, new_vcf_record.REF, new_vcf_record.ALT)
        expected = _MockVcfRecord(POS = 5, REF = 'C', ALT = ['G'])

        self.assertEqual(expected, result)


    def test_StartsAtNonSite_EndsAtSite(self):
        # base sequence:      T TAT CGG
        # secondary sequence: T G   CGG
        secondary_reference_length = 5
        base_records = [
            _MockVcfRecord(POS=2, REF="TAT", ALT=['G'])
        ]
        secondary_regions = discover._flag_personalised_reference_regions(base_records, secondary_reference_length)

        secondary_vcf_record = _MockVcfRecord(POS=1, REF='TG', ALT=['TAA'])
        new_vcf_record = discover._rebase_vcf_record(secondary_vcf_record, secondary_regions)

        result = _MockVcfRecord(new_vcf_record.POS, new_vcf_record.REF, new_vcf_record.ALT)
        expected = _MockVcfRecord(1, 'TTAT', ['TAA'])

        self.assertEqual(expected, result)


    def test_SiteInBetweenNonSites(self):
        """
        A test case where the variation on top of the inferred reference overlaps: a non-variant site, a variant site,
        and a non-variant site in the prg.

        What we need is for the rebased REF to include all three sites.
        """
        # base sequ: T TAT CGG
        # secondary: T G   CGG
        secondary_reference_length = 5
        base_records = [
            _MockVcfRecord(POS=2, REF="TAT", ALT=['G'])
        ]

        secondary_vcf_record = _MockVcfRecord(POS=1, REF='TGCG', ALT=['GGCT'])

        secondary_regions = discover._flag_personalised_reference_regions(base_records, secondary_reference_length)
        new_vcf_record = discover._rebase_vcf_record(secondary_vcf_record, secondary_regions)

        result = _MockVcfRecord(new_vcf_record.POS, new_vcf_record.REF, new_vcf_record.ALT)
        expected = _MockVcfRecord(POS = 1, REF = 'TTATCG', ALT = ['GGCT'])

        self.assertEqual(expected, result)


    def test_SNP_OnTopOfIndel(self):
        """
        A test case where we find a SNP on top of an insertion in the inferred reference.

        What we need is for the rebased ALT to include the flanking ALT bases, which are implied to be present in the secondary_vcf_record.
        """
        # base sequ: T TAT CGG T     A
        # secondary: T G   CGG TCTGC A
        secondary_reference_length = 8
        base_records = [
            _MockVcfRecord(POS=2, REF="TAT", ALT=['G']),
            _MockVcfRecord(POS=8, REF="T", ALT=['TCTGC'])
        ]

        secondary_vcf_record = _MockVcfRecord(POS=9, REF='G', ALT=['A'])

        secondary_regions = discover._flag_personalised_reference_regions(base_records, secondary_reference_length)
        new_vcf_record = discover._rebase_vcf_record(secondary_vcf_record, secondary_regions)

        result = _MockVcfRecord(new_vcf_record.POS, new_vcf_record.REF, new_vcf_record.ALT)
        expected = _MockVcfRecord(8,'T',['TCTAC'])

        self.assertEqual(expected, result)

    def test_Deletion_OnTopOfDeletion_WithExtraDeletionInNonSite(self):
        """
        A test case where we discover a deletion on top of a deletion in a variant site;
        as well as an extra deletion in a non-variant site.

        There is also a SNP among the original deletion, to make it plausible that quasimap/infer picks this variant.

        To make it harder, the discovered variation is also reported inside a variant site, so we expect the rebased ALT to be elongated.

        We expect the rebased REF to include all deleted bases.
        """
        # base reference:     CAA C GCTA CAA
        # inferred reference: C   C GAT  CAA


        secondary_reference_length = 8
        base_records = [
           _MockVcfRecord(POS=1, REF="CAA", ALT=['C']),
           _MockVcfRecord(POS=5, REF="GCTA", ALT=['GAT'])
        ]

        secondary_vcf_record = _MockVcfRecord(POS = 4, REF='ATC', ALT=['A'])


        secondary_regions = discover._flag_personalised_reference_regions(base_records, secondary_reference_length)
        new_vcf_record = discover._rebase_vcf_record(secondary_vcf_record, secondary_regions)

        result = _MockVcfRecord(new_vcf_record.POS, new_vcf_record.REF, new_vcf_record.ALT)
        expected = _MockVcfRecord(POS = 5, REF = 'GCTAC', ALT = ['GA'])

        self.assertEqual(expected, result)



class TestGetBaseReference(unittest.TestCase):
    """
    Pick the first allele of each variant site of the prg to make 'base' reference.
    """
    def test_GivenPrgWithTwoSites_CorrectReferenceExtracted(self):
        test_prg = "test_prg"
        test_inferred = "test_inferred"

        prg_seq = 'TT5A6T5AA7C8A7'

        allele_indexes = (0 for _ in range(4))

        with open(test_prg, "w") as prg:
            prg.write(prg_seq)


        prg_parser = prg_local_parser.Prg_Local_Parser(test_prg, test_inferred, fasta_header= '', allele_indexes = allele_indexes)
        prg_parser.parse()

        with open(test_inferred,"r") as inferred:
            result = inferred.readlines()[-1].strip()

        expected = 'TTAAAC'
        os.remove(test_prg); os.remove(test_inferred)

        self.assertEqual(expected, result)


class _Sample:
    def __init__(self,genotype):
        self.gt_alleles = genotype


class _MockVcfRecord:
    def __init__(self, POS, REF, ALT, samples=[]):
        self.POS = POS
        self.REF = REF
        self.ALT = [vcf.model._Substitution(x) for x in ALT]
        self.samples = samples

    def __repr__(self):
        return str(self.__dict__)

    def __eq__(self, other):
        return self.POS == other.POS \
                and self.REF == other.REF \
                and self.ALT == other.ALT


if __name__ == "__main__":
    unittest.main()