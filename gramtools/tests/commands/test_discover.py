import unittest
import vcf
import os

from ...commands import discover
from ... import prg_local_parser


class TestSecondaryRegionsForBaseSites(unittest.TestCase):
    def test_SingleBaseAlt_CorrectRegion(self):
        # base sequence:      T TAT CGG
        # secondary sequence: T G   CGG
        base_records = [
            _MockVcfRecord(POS=2, REF="TAT", ALT=['G'])
        ]
        result = discover.mark_base_site_regions(base_records)

        expected = [
            discover._Region(start=1, end=1, record=base_records[0], offset=None)
        ]
        self.assertEqual(expected, result)

    def test_AltLongerThanRef_CorrectRegion(self):
        # base sequence:      T TAT    CGG
        # secondary sequence: T GCCAC  CGG
        base_records = [
            _MockVcfRecord(POS=2, REF="TAT", ALT=['GCCAC'])
        ]
        result = discover.mark_base_site_regions(base_records)

        expected = [
            discover._Region(start=1, end=5, record=base_records[0], offset=None)
        ]
        self.assertEqual(expected, result)

    def test_TwoRecords_CorrectRegions(self):
        # base sequence:      T TAT    C G   G
        # secondary sequence: T GCCAC  C TTT G
        base_records = [
            _MockVcfRecord(POS=2, REF="TAT", ALT=['GCCAC']),
            _MockVcfRecord(POS=6, REF="G", ALT=['TTT'])
        ]
        result = discover.mark_base_site_regions(base_records)

        expected = [
            discover._Region(start=1, end=5, record=base_records[0], offset=None),
            discover._Region(start=7, end=9, record=base_records[1], offset=None),
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
        result = discover.mark_base_site_regions(base_records)

        expected = [
            discover._Region(start=1, end=5, record=base_records[0], offset=None),
            discover._Region(start=6, end=8, record=base_records[1], offset=None),
            discover._Region(start=9, end=10, record=base_records[2], offset=None),
        ]
        self.assertEqual(expected, result)


class TestSecondaryRegionsForBaseNonsites(unittest.TestCase):
    def test_SingleBaseAlt_CorrectRegions(self):
        # base sequence:      T TAT CGG
        # secondary sequence: T G   CGG
        secondary_reference_length = 5
        base_records = [
            _MockVcfRecord(POS=2, REF="TAT", ALT=['G'])
        ]
        secondary_site_regions = discover.mark_base_site_regions(base_records)
        result = discover.mark_base_nonsite_regions(secondary_site_regions,
                                                    base_records,
                                                    secondary_reference_length)
        expected = [
            discover._Region(start=0, end=0, record=None, offset=0),
            discover._Region(start=2, end=4, record=None, offset=2),
        ]
        self.assertEqual(expected, result)

    def test_AltLongerThanRef_CorrectRegion(self):
        # base sequence:      T TAT    CGG
        # secondary sequence: T GCCAC  CGG
        secondary_reference_length = 9
        base_records = [
            _MockVcfRecord(POS=2, REF="TAT", ALT=['GCCAC'])
        ]
        secondary_site_regions = discover.mark_base_site_regions(base_records)
        result = discover.mark_base_nonsite_regions(secondary_site_regions,
                                                    base_records,
                                                    secondary_reference_length)
        expected = [
            discover._Region(start=0, end=0, record=None, offset=0),
            discover._Region(start=6, end=8, record=None, offset=-2),
        ]
        self.assertEqual(expected, result)

    def test_TwoRecords_CorrectRegions(self):
        # base sequence:      T TAT    C G   G
        # secondary sequence: T GCCAC  C TTT G
        secondary_reference_length = 11
        base_records = [
            _MockVcfRecord(POS=2, REF="TAT", ALT=['GCCAC']),
            _MockVcfRecord(POS=6, REF="G", ALT=['TTT'])
        ]
        secondary_site_regions = discover.mark_base_site_regions(base_records)
        result = discover.mark_base_nonsite_regions(secondary_site_regions,
                                                    base_records,
                                                    secondary_reference_length)
        expected = [
            discover._Region(start=0, end=0, record=None, offset=0),
            discover._Region(start=6, end=6, record=None, offset=-2),
            discover._Region(start=10, end=10, record=None, offset=-4),
        ]
        self.assertEqual(expected, result)

    def test_ThreeAdjacentRecords_CorrectRegions(self):
        # base sequence:      T TAT    C   G  G
        # secondary sequence: T GCCAC  TCT AA G
        secondary_reference_length = 12
        base_records = [
            _MockVcfRecord(POS=2, REF="TAT", ALT=['GCCAC']),
            _MockVcfRecord(POS=5, REF="C", ALT=['TCT']),
            _MockVcfRecord(POS=6, REF="G", ALT=['AA']),
        ]
        secondary_site_regions = discover.mark_base_site_regions(base_records)
        result = discover.mark_base_nonsite_regions(secondary_site_regions,
                                                    base_records,
                                                    secondary_reference_length)
        expected = [
            discover._Region(start=0, end=0, record=None, offset=0),
            discover._Region(start=11, end=11, record=None, offset=-5),
        ]
        self.assertEqual(expected, result)

    def test_NoRecords_NoRegions(self):
        # base sequence:      TTATCGG
        # secondary sequence: TTATCGG
        secondary_reference_length = 7
        base_records = []
        secondary_site_regions = discover.mark_base_site_regions(base_records)
        result = discover.mark_base_nonsite_regions(secondary_site_regions,
                                                    base_records,
                                                    secondary_reference_length)
        expected = []
        self.assertEqual(expected, result)

    def test_OnlyAdjacentRecords_NoNonsiteRegions(self):
        # base sequence:      T TAT    C   G  G
        # secondary sequence: A GCCAC  TCT AA T
        secondary_reference_length = 12
        base_records = [
            _MockVcfRecord(POS=1, REF="T", ALT=['A']),
            _MockVcfRecord(POS=2, REF="TAT", ALT=['GCCAC']),
            _MockVcfRecord(POS=5, REF="C", ALT=['TCT']),
            _MockVcfRecord(POS=6, REF="G", ALT=['AA']),
            _MockVcfRecord(POS=7, REF="G", ALT=['T']),
        ]
        secondary_site_regions = discover.mark_base_site_regions(base_records)
        result = discover.mark_base_nonsite_regions(secondary_site_regions,
                                                    base_records,
                                                    secondary_reference_length)
        expected = []
        self.assertEqual(expected, result)

    def test_EndsInRecord_CorrectRegions(self):
        # base sequence:      TTATC GG
        # secondary sequence: TTATC TAC
        secondary_reference_length = 8
        base_records = [
            _MockVcfRecord(POS=6, REF="GG", ALT=['TAC']),
        ]
        secondary_site_regions = discover.mark_base_site_regions(base_records)
        result = discover.mark_base_nonsite_regions(secondary_site_regions,
                                                    base_records,
                                                    secondary_reference_length)
        expected = [
            discover._Region(start=0, end=4, record=None, offset=0)
        ]
        self.assertEqual(expected, result)

    def test_StartsInRecord_CorrectRegions(self):
        # base sequence:      TTA TCGG
        # secondary sequence: TC  TCGG
        secondary_reference_length = 6
        base_records = [
            _MockVcfRecord(POS=1, REF="TTA", ALT=['TC']),
        ]
        secondary_site_regions = discover.mark_base_site_regions(base_records)
        result = discover.mark_base_nonsite_regions(secondary_site_regions,
                                                    base_records,
                                                    secondary_reference_length)
        expected = [
            discover._Region(start=2, end=5, record=None, offset=1)
        ]
        self.assertEqual(expected, result)


class TestIterPairs(unittest.TestCase):
    def test_EvenNumberOfElements_CorrectPairing(self):
        elements = [1, 2, 3, 4]
        result = list(discover.IterPairs(elements))

        expected = [
            (None, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, None),
        ]
        self.assertEqual(expected, result)

    def test_OddNumberOfElements_CorrectPairing(self):
        elements = [1, 2, 3, 4, 5]
        result = list(discover.IterPairs(elements))

        expected = [
            (None, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, None),
        ]
        self.assertEqual(expected, result)

    def test_NoElements_NoPairs(self):
        elements = []
        result = list(discover.IterPairs(elements))

        expected = []
        self.assertEqual(expected, result)

    def test_OneElement_CorrectPair(self):
        elements = [1]
        result = list(discover.IterPairs(elements))

        expected = [
            (None, 1),
            (1, None),
        ]
        self.assertEqual(expected, result)

    def test_TwoElements_CorrectPairs(self):
        elements = [1, 2]
        result = list(discover.IterPairs(elements))

        expected = [
            (None, 1),
            (1, 2),
            (2, None),
        ]
        self.assertEqual(expected, result)


class TestMergeRegions(unittest.TestCase):
    def test_SingleBaseRegionMiddleMerge_CorrectMerge(self):
        site_regions = [
            discover._Region(start=2, end=5, record=None, offset=None),
            discover._Region(start=7, end=7, record=None, offset=None),
        ]
        nonsite_regions = [
            discover._Region(start=6, end=6, record=None, offset=None),
        ]
        result = discover._merge_regions(site_regions, nonsite_regions)

        expected = [
            discover._Region(start=2, end=5, record=None, offset=None),
            discover._Region(start=6, end=6, record=None, offset=None),
            discover._Region(start=7, end=7, record=None, offset=None),
        ]
        self.assertEqual(expected, result)

    def test_AlternatingRegions_CorrectMerge(self):
        site_regions = [
            discover._Region(start=2, end=5, record=None, offset=None),
            discover._Region(start=7, end=7, record=None, offset=None),
        ]
        nonsite_regions = [
            discover._Region(start=6, end=6, record=None, offset=None),
            discover._Region(start=8, end=10, record=None, offset=None),
        ]
        result = discover._merge_regions(site_regions, nonsite_regions)

        expected = [
            discover._Region(start=2, end=5, record=None, offset=None),
            discover._Region(start=6, end=6, record=None, offset=None),
            discover._Region(start=7, end=7, record=None, offset=None),
            discover._Region(start=8, end=10, record=None, offset=None),
        ]
        self.assertEqual(expected, result)


class TestSearchRegions(unittest.TestCase):
    def test_RecordStartsInRangeMid_CorrectRegionIndex(self):
        secondary_regions = [
            discover._Region(start=2, end=5, record=None, offset=None),
            discover._Region(start=6, end=6, record=None, offset=None),
            discover._Region(start=7, end=10, record=None, offset=None),
            discover._Region(start=11, end=12, record=None, offset=None),
        ]
        vcf_record = _MockVcfRecord(POS=10, REF='A', ALT=['A'])
        result = discover._find_start_region_index(vcf_record, secondary_regions)

        expected = 2
        self.assertEqual(expected, result)

    def test_RecordStartsAtFirstRegion_CorrectRegionIndex(self):
        secondary_regions = [
            discover._Region(start=2, end=5, record=None, offset=None),
            discover._Region(start=6, end=6, record=None, offset=None),
            discover._Region(start=7, end=10, record=None, offset=None),
            discover._Region(start=11, end=12, record=None, offset=None),
        ]
        vcf_record = _MockVcfRecord(POS=3, REF='A', ALT=['A'])
        result = discover._find_start_region_index(vcf_record, secondary_regions)

        expected = 0
        self.assertEqual(expected, result)

    def test_RecordStartsAtSecondRegion_CorrectRegionIndex(self):
        secondary_regions = [
            discover._Region(start=2, end=5, record=None, offset=None),
            discover._Region(start=6, end=6, record=None, offset=None),
            discover._Region(start=7, end=10, record=None, offset=None),
            discover._Region(start=11, end=12, record=None, offset=None),
        ]
        vcf_record = _MockVcfRecord(POS=7, REF='A', ALT=['A'])
        result = discover._find_start_region_index(vcf_record, secondary_regions)

        expected = 1
        self.assertEqual(expected, result)


class TestFindOverlapRegions(unittest.TestCase):
    def test_RecordEncapsulatedInFirstRegion_ReturnFirstRegion(self):
        secondary_regions = [
            discover._Region(start=2, end=6, record=None, offset=None),
            discover._Region(start=7, end=7, record=None, offset=None),
        ]
        vcf_record = _MockVcfRecord(POS=4, REF='AA', ALT=['A'])
        first_region_index = discover._find_start_region_index(vcf_record, secondary_regions)
        result = discover._find_overlap_regions(vcf_record, first_region_index, secondary_regions)

        expected = [
            discover._Region(start=2, end=6, record=None, offset=None),
        ]
        self.assertEqual(expected, result)

    def test_RecordReachFirstRegionBoundary_ReturnFirstRegion(self):
        secondary_regions = [
            discover._Region(start=2, end=5, record=None, offset=None),
            discover._Region(start=6, end=6, record=None, offset=None),
        ]
        vcf_record = _MockVcfRecord(POS=3, REF='AAAA', ALT=['A'])
        first_region_index = discover._find_start_region_index(vcf_record, secondary_regions)
        result = discover._find_overlap_regions(vcf_record, first_region_index, secondary_regions)

        expected = [
            discover._Region(start=2, end=5, record=None, offset=None),
        ]
        self.assertEqual(expected, result)

    def test_RecordBarelyOverlapSecondRegion_CorrectRegions(self):
        secondary_regions = [
            discover._Region(start=2, end=5, record=None, offset=None),
            discover._Region(start=6, end=6, record=None, offset=None),
            discover._Region(start=7, end=10, record=None, offset=None),
        ]
        vcf_record = _MockVcfRecord(POS=3, REF='AAAAC', ALT=['A'])
        first_region_index = discover._find_start_region_index(vcf_record, secondary_regions)
        result = discover._find_overlap_regions(vcf_record, first_region_index, secondary_regions)

        expected = [
            discover._Region(start=2, end=5, record=None, offset=None),
            discover._Region(start=6, end=6, record=None, offset=None),
        ]
        self.assertEqual(expected, result)

    def test_RecordEncapsulatedInSecondRegion_ReturnSecondRegion(self):
        secondary_regions = [
            discover._Region(start=2, end=5, record=None, offset=None),
            discover._Region(start=6, end=6, record=None, offset=None),
            discover._Region(start=7, end=10, record=None, offset=None),
        ]
        vcf_record = _MockVcfRecord(POS=7, REF='A', ALT=['A'])
        first_region_index = discover._find_start_region_index(vcf_record, secondary_regions)
        result = discover._find_overlap_regions(vcf_record, first_region_index, secondary_regions)

        expected = [
            discover._Region(start=6, end=6, record=None, offset=None),
        ]
        self.assertEqual(expected, result)

    def test_RecordOverlapLastRegion_ReturnLastRegion(self):
        secondary_regions = [
            discover._Region(start=2, end=5, record=None, offset=None),
            discover._Region(start=6, end=6, record=None, offset=None),
            discover._Region(start=7, end=10, record=None, offset=None),
        ]
        vcf_record = _MockVcfRecord(POS=10, REF='AA', ALT=['A'])
        first_region_index = discover._find_start_region_index(vcf_record, secondary_regions)
        result = discover._find_overlap_regions(vcf_record, first_region_index, secondary_regions)

        expected = [
            discover._Region(start=7, end=10, record=None, offset=None),
        ]
        self.assertEqual(expected, result)


class TestRebaseVcfRecord(unittest.TestCase):
    def test_NoSiteOverlap_CorrectPos(self):
        # base sequence:      T TAT CGG
        # secondary sequence: T G   CGG
        secondary_reference_length = 5
        base_records = [
            _MockVcfRecord(POS=2, REF="TAT", ALT=['G'])
        ]
        secondary_regions = discover._flag_personalised_reference_regions(base_records, secondary_reference_length)
        secondary_vcf_record = _MockVcfRecord(POS=3, REF='C', ALT=['G'])
        new_vcf_record = discover._rebase_vcf_record(secondary_vcf_record, secondary_regions)
        result = new_vcf_record.POS

        expected = 5
        self.assertEqual(expected, result)

    def test_StartsAtNonSite_PosUnmodified(self):
        # base sequence:      T TAT CGG
        # secondary sequence: T G   CGG
        secondary_reference_length = 5
        base_records = [
            _MockVcfRecord(POS=2, REF="TAT", ALT=['G'])
        ]
        secondary_regions = discover._flag_personalised_reference_regions(base_records, secondary_reference_length)
        secondary_vcf_record = _MockVcfRecord(POS=1, REF='TG', ALT=['TAA'])
        new_vcf_record = discover._rebase_vcf_record(secondary_vcf_record, secondary_regions)
        result = new_vcf_record.POS

        expected = 1
        self.assertEqual(expected, result)


class TestPadVcfRecordStart(unittest.TestCase):
    def test_RecordStartsOneBaseIntoSite_CorrectModifiedRecord(self):
        # base sequence:      T TAT  CGG
        # secondary sequence: T GAAT CGG
        base_records = [
            _MockVcfRecord(POS=2, REF="TAT", ALT=['GAAT'])
        ]
        regions = [
            discover._Region(start=1, end=4, record=base_records[0], offset=None)
        ]

        vcf_record = _MockVcfRecord(POS=3, REF="AA", ALT=['C'])
        result = discover._pad_vcf_record_start(vcf_record, regions)

        expected = _MockVcfRecord(POS=2, REF='TAT', ALT=['GC'])
        self.assertEqual(expected, result)

    def test_RecordStartsAtSiteEndBase_CorrectModifiedRecord(self):
        # base sequence:      T TAT  CGG
        # secondary sequence: T GAAT CGG
        base_records = [
            _MockVcfRecord(POS=2, REF="TAT", ALT=['GAAT'])
        ]
        regions = [
            discover._Region(start=1, end=4, record=base_records[0], offset=None)
        ]

        vcf_record = _MockVcfRecord(POS=4, REF="AT", ALT=['G'])
        result = discover._pad_vcf_record_start(vcf_record, regions)

        expected = _MockVcfRecord(POS=2, REF='TAT', ALT=['GAG'])
        self.assertEqual(expected, result)

    def test_RecordStartsAtSiteStart_CorrectModifiedRecord(self):
        # base sequence:      T TAT  CGG
        # secondary sequence: T GAAT CGG
        base_records = [
            _MockVcfRecord(POS=2, REF="TAT", ALT=['GAAT'])
        ]
        regions = [
            discover._Region(start=1, end=4, record=base_records[0], offset=None)
        ]

        vcf_record = _MockVcfRecord(POS=2, REF="GA", ALT=['C'])
        result = discover._pad_vcf_record_start(vcf_record, regions)

        expected = _MockVcfRecord(POS=2, REF='TAT', ALT=['C'])
        self.assertEqual(expected, result)


class TestPadVcfRecordEnd(unittest.TestCase):
    def test_RecordStartsWithinSite_AltExtendedWithSiteEnd(self):
        # base sequence:      T TAT  CGG
        # secondary sequence: T GAAT CGG
        base_records = [
            _MockVcfRecord(POS=2, REF="TAT", ALT=['GAAT'])
        ]
        regions = [
            discover._Region(start=1, end=4, record=base_records[0], offset=None)
        ]

        vcf_record = _MockVcfRecord(POS=3, REF="AA", ALT=['C'])
        result = discover._pad_vcf_record_end(vcf_record, regions)

        expected = _MockVcfRecord(POS=3, REF="AA", ALT=['CT'])
        self.assertEqual(expected, result)

    def test_RecordStartsBeforeSite_AltExtendedWithSiteEnd(self):
        # base sequence:      T TAT     CGG
        # secondary sequence: T GAATTTG CGG
        base_records = [
            _MockVcfRecord(POS=2, REF="TAT", ALT=['GAATTTG'])
        ]
        regions = [
            discover._Region(start=0, end=0, record=None, offset=0),
            discover._Region(start=1, end=7, record=base_records[0], offset=None)
        ]

        vcf_record = _MockVcfRecord(POS=1, REF="TGAA", ALT=['CAGT'])
        result = discover._pad_vcf_record_end(vcf_record, regions)

        expected = _MockVcfRecord(POS=1, REF="TGAA", ALT=['CAGTTTTG'])
        self.assertEqual(expected, result)

    def test_AltShorterThanRef_AltExtendedWithSiteEnd(self):
        # base sequence:      T TAT     CGG
        # secondary sequence: T GAATTTG CGG
        base_records = [
            _MockVcfRecord(POS=2, REF="TAT", ALT=['GAATTTG'])
        ]
        regions = [
            discover._Region(start=0, end=0, record=None, offset=0),
            discover._Region(start=1, end=7, record=base_records[0], offset=None)
        ]

        vcf_record = _MockVcfRecord(POS=1, REF="TGAA", ALT=['C'])
        result = discover._pad_vcf_record_end(vcf_record, regions)

        expected = _MockVcfRecord(POS=1, REF="TGAA", ALT=['CTTTG'])
        self.assertEqual(expected, result)

    def test_AltLongerThanRef_AltExtendedWithSiteEnd(self):
        # base sequence:      TACTAC TAT  CGG
        # secondary sequence: TACTAC GAAT CGG
        base_records = [
            _MockVcfRecord(POS=7, REF="TAT", ALT=['GAAT'])
        ]
        regions = [
            discover._Region(start=0, end=5, record=None, offset=0),
            discover._Region(start=6, end=9, record=base_records[0], offset=None)
        ]

        vcf_record = _MockVcfRecord(POS=4, REF="TACG", ALT=['CCCTTT'])
        result = discover._pad_vcf_record_end(vcf_record, regions)
        print(result)

        expected = _MockVcfRecord(POS=4, REF="TACG", ALT=['CCCTTTAAT'])
        self.assertEqual(expected, result)


class TestGetReference(unittest.TestCase):
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


class TestGetInferredReference(unittest.TestCase):

    def setUp(self):
        self.test_ref = "test_ref.fa"
        self.test_inferred = "test_inferred.fa"

    def tearDown(self):
        os.remove(self.test_ref)
        os.remove(self.test_inferred)


    def test_VcfEntryAltShorterThanRef_CorrectInferReference(self):

        with open(self.test_ref, "w") as ref:
            ref.write('>First_line \n TTATCGG')

        vcf_reader = iter([
            _MockVcfRecord(POS=2, REF="TAT", ALT=['G'], samples = [_Sample(genotype = [1])])
        ])
        discover._produce_inferred_reference(self.test_ref, vcf_reader, self.test_inferred, description='')

        with open(self.test_inferred, "r") as inferred:
            result = inferred.readlines()[-1].strip()

        expected = 'TGCGG'
        self.assertEqual(expected, result)

    def test_VcfEntryAltLongerThanRef_CorrectInferReference(self):

        with open(self.test_ref, "w") as ref:
            ref.write('>First_line \n TTATCGG')

        vcf_reader = iter([
            _MockVcfRecord(POS=2, REF="TAT", ALT=['GCG'], samples = [_Sample(genotype = [1])]),
        ])
        discover._produce_inferred_reference(self.test_ref, vcf_reader, self.test_inferred, description = '')

        with open(self.test_inferred, "r") as inferred:
            result = inferred.readlines()[-1].strip()

        expected = 'TGCGCGG'
        self.assertEqual(expected, result)

    def test_MultipleVcfEntries_CorrectInferReference(self):

        with open(self.test_ref, "w") as ref:
            ref.write('>First_line \n TTATCGG')

        vcf_reader = iter([
            _MockVcfRecord(POS=2, REF="TAT", ALT=['GCG'], samples = [_Sample(genotype = [1])]),
            _MockVcfRecord(POS=6, REF="G", ALT=['A'], samples = [_Sample(genotype = [1])]),
        ])
        discover._produce_inferred_reference(self.test_ref, vcf_reader, self.test_inferred, description = '')

        with open(self.test_inferred, "r") as inferred:
            result = inferred.readlines()[-1].strip()

        expected = 'TGCGCAG'
        self.assertEqual(expected, result)

    def test_VcfEntryModifiesLastBase_CorrectInferReference(self):

        with open(self.test_ref, "w") as ref:
            ref.write('>First_line \n TTATCGG')

        vcf_reader = iter([
            _MockVcfRecord(POS=7, REF="G", ALT=['GCG'],samples = [_Sample(genotype = [1])]),
        ])
        discover._produce_inferred_reference(self.test_ref, vcf_reader, self.test_inferred, description='')

        with open(self.test_inferred, "r") as inferred:
            result = inferred.readlines()[-1].strip()

        expected = 'TTATCGGCG'

        self.assertEqual(expected, result)

if __name__ == "__main__":
    unittest.main()