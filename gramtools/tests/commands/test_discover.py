import unittest
import collections

from ...commands import discover


class TestGetReference(unittest.TestCase):
    def test_GivenPrgWithTwoSites_CorrectReferenceExtracted(self):
        prg_seq = 'TT5A6T5AA7C8A7'
        result = discover._get_reference(prg_seq)

        expected = ['T', 'T', 'A', 'A', 'A', 'C']
        self.assertEqual(expected, result)


_MockVcfRecord = collections.namedtuple('MockVcfRecord', ['POS', 'REF', 'ALT'])


class TestGetInferredReference(unittest.TestCase):
    def test_VcfEntryAltShorterThanRef_CorrectInferReference(self):
        reference = ['T', 'T', 'A', 'T', 'C', 'G', 'G']
        vcf_reader = iter([
            _MockVcfRecord(POS=2, REF="TAT", ALT=['G'])
        ])
        result = discover._get_inferred_reference(reference, vcf_reader)

        expected = ['T', 'G', 'C', 'G', 'G']
        self.assertEqual(expected, result)

    def test_VcfEntryAltLongerThanRef_CorrectInferReference(self):
        reference = ['T', 'T', 'A', 'T', 'C', 'G', 'G']
        vcf_reader = iter([
            _MockVcfRecord(POS=2, REF="TAT", ALT=['GCG']),
        ])
        result = discover._get_inferred_reference(reference, vcf_reader)

        expected = ['T', 'G', 'C', 'G', 'C', 'G', 'G']
        self.assertEqual(expected, result)

    def test_MultipleVcfEntries_CorrectInferReference(self):
        reference = ['T', 'T', 'A', 'T', 'C', 'G', 'G']
        vcf_reader = iter([
            _MockVcfRecord(POS=2, REF="TAT", ALT=['GCG']),
            _MockVcfRecord(POS=6, REF="G", ALT=['A']),
        ])
        result = discover._get_inferred_reference(reference, vcf_reader)

        expected = ['T', 'G', 'C', 'G', 'C', 'A', 'G']
        self.assertEqual(expected, result)

    def test_VcfEntryModifiesLastBase_CorrectInferReference(self):
        reference = ['T', 'T', 'A', 'T', 'C', 'G', 'G']
        vcf_reader = iter([
            _MockVcfRecord(POS=7, REF="G", ALT=['GCG']),
        ])
        result = discover._get_inferred_reference(reference, vcf_reader)

        expected = ['T', 'T', 'A', 'T', 'C', 'G', 'G', 'C', 'G']
        self.assertEqual(expected, result)
