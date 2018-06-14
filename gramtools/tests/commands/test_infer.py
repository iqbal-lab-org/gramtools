import copy
import unittest
from unittest.mock import Mock

from ...commands import infer


class TestIsInt(unittest.TestCase):
    def test_MultipleCases_CorrectResults(self):
        cases = {
            'A': False,
            ' ': False,
            '': False,
            None: False,

            '123': True,
            '42': True,
        }
        for value, expected in cases.items():
            result = infer._is_int(value)
            self.assertEqual(result, expected)


class TestParsePrgChars(unittest.TestCase):
    def test_MultidigitNumbers_Correct(self):
        chars = 'GT123ACT321'
        result = list(infer._parse_prg_chars(chars))
        expected = ['G', 'T', '123', 'A', 'C', 'T', '321']
        self.assertEqual(result, expected)

    def test_MultidigitNumber_Correct(self):
        chars = '123'
        result = list(infer._parse_prg_chars(chars))
        expected = ['123']
        self.assertEqual(result, expected)

    def test_MultidigitNumbersAtStartAndEnd_Correct(self):
        chars = '123T321A42'
        result = list(infer._parse_prg_chars(chars))
        expected = ['123', 'T', '321', 'A', '42']
        self.assertEqual(result, expected)


class TestParsePrgStructure(unittest.TestCase):
    def test_GivenPrg_CorrectCharsInCursor(self):
        chars = 'ACGT5AA6T5C'
        cursor = infer._Cursor()
        result = []
        for cursor in infer._parse_prg_structure(chars, cursor):
            result.append(cursor.char)
        expected = ['A', 'C', 'G', 'T', '5', 'A', 'A', '6', 'T', '5', 'C']
        self.assertEqual(result, expected)

    def test_GivenPrg_CorrectJustLeftSiteFlagInCursor(self):
        chars = 'ACGT5AA6T5C'
        cursor = infer._Cursor()
        result = []
        for cursor in infer._parse_prg_structure(chars, cursor):
            result.append(cursor.just_left_site)
        expected = [False, False, False, False, False, False, False, False, False, True, False]
        self.assertEqual(result, expected)

    def test_GivenPrg_CorrectOnMarkerFlagInCursor(self):
        chars = 'ACGT5AA6T5C'
        cursor = infer._Cursor()
        result = []
        for cursor in infer._parse_prg_structure(chars, cursor):
            result.append(cursor.on_marker)
        expected = [False, False, False, False, True, False, False, True, False, True, False]
        self.assertEqual(result, expected)

    def test_GivenPrgOneSite_CorrectAlleleIdInCursor(self):
        chars = 'ACGT5AA6T5C'
        cursor = infer._Cursor()
        result = []
        for cursor in infer._parse_prg_structure(chars, cursor):
            result.append(cursor.allele_id)
        expected = [
            None, None, None, None, # ACGT
            None,                   # 5
            0, 0,                   # AA
            None,                   # 6
            1,                      # T
            None, None              # 5C
        ]
        self.assertEqual(result, expected)

    def test_GivenPrgTwoSites_CorrectAlleleIdInCursor(self):
        chars = 'ACGT5AA6T5C7AA8T7'
        cursor = infer._Cursor()
        result = []
        for cursor in infer._parse_prg_structure(chars, cursor):
            result.append(cursor.allele_id)
        expected = [
            None, None, None, None, # ACGT
            None,                   # 5
            0, 0,                   # AA
            None,                   # 6
            1,                      # T
            None, None,             # 5C
            None,                   # 7
            0, 0,                   # AA
            None,                   # 8
            1,                      # T
            None,                   # 7
        ]
        self.assertEqual(expected, result)

    def test_GivenPrg_CorrectSiteMarkerInCursor(self):
        chars = 'ACGT5AA6T5C'
        cursor = infer._Cursor()
        result = []
        for cursor in infer._parse_prg_structure(chars, cursor):
            result.append(cursor.site_marker)
        expected = [
            None, None, None, None, # ACGT
            5,                      # 5
            5, 5,                   # AA
            5,                      # 6
            5,                      # T
            5,                      # 5
            None,                   # C
        ]
        self.assertEqual(result, expected)


class TestReadChunk(unittest.TestCase):
    def test_PrgEndsInNumber_CorrectChunk(self):
        file_handle = Mock()
        file_handle.read.side_effect = ['ACGT5', '']
        result = infer._read_chunk(file_handle)
        expected = 'ACGT5'
        self.assertEqual(expected, result)

    def test_FirstChunkEndsInNumber_ConsumesCharectersUntilNonNumber(self):
        file_handle = Mock()
        file_handle.read.side_effect = ['ACGT5', '6', 'C', 'A']
        result = infer._read_chunk(file_handle)
        expected = 'ACGT56C'
        self.assertEqual(expected, result)

    def test_PrgEndsInNumber_ReadCalledTwice(self):
        file_handle = Mock()
        file_handle.read.side_effect = ['ACGT5', '']
        infer._read_chunk(file_handle)

        result = file_handle.read.call_count
        expected = 2
        self.assertEqual(expected, result)

    def test_PrgNotEndInNumber_ReadCalledOnce(self):
        file_handle = Mock()
        file_handle.read.side_effect = ['ACGT5C']
        infer._read_chunk(file_handle)

        result = file_handle.read.call_count
        expected = 1
        self.assertEqual(expected, result)


class TestAssembleReference(unittest.TestCase):
    def test_PrgEndsWithNumber_CorrectReference(self):
        chars = 'AC5G6AA5CCC7AC8T7'
        cursor = infer._Cursor()
        prg_parser = [copy.copy(cursor)
                      for cursor in infer._parse_prg_structure(chars, cursor)]
        allele_indexes = iter([1, 0])
        cache_writer = []
        infer._dump_fasta(prg_parser, allele_indexes, cache_writer)

        result = ''.join(cache_writer)
        expected = 'ACAACCCAC'
        self.assertEqual(expected, result)

    def test_GivenPrgAndAlleleIndexes_CorrectReference(self):
        chars = 'AC5G6AA5CCC7AC8T7TT'
        cursor = infer._Cursor()
        prg_parser = [copy.copy(cursor)
                      for cursor in infer._parse_prg_structure(chars, cursor)]
        allele_indexes = iter([1, 0])
        cache_writer = []
        infer._dump_fasta(prg_parser, allele_indexes, cache_writer)

        result = ''.join(cache_writer)
        expected = 'ACAACCCACTT'
        self.assertEqual(expected, result)
