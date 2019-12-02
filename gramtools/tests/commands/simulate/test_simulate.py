import unittest

from gramtools.commands.simulate import prg_regions_parser
from gramtools.commands.simulate import simulate
from ... import common


class TestGenerateReads(unittest.TestCase):
    def _analyse_case(self, read_length, prg_structure,
                      expected, max_num_reads=None):
        prg_seq = common.compose_prg(prg_structure)
        regions = prg_regions_parser.parse(prg_seq)
        reads = set(simulate._generate_reads(read_length, regions,
                                             max_num_reads))
        self.assertEqual(reads, expected)

    def test_variantSiteManyLengths_readLengthsCorrect(self):
        prg_structure = [
            ['AGAGAG'],
            ['TC', 'CATT'],
        ]
        expected = {'GAGTC', 'GCATT'}
        read_length = 5
        self._analyse_case(read_length, prg_structure, expected)

    def test_minimalReadLengths_correctReads(self):
        prg_structure = [
            ['AGCGAG'],
            ['GC', 'AT'],
        ]
        expected = {'C', 'T'}
        read_length = 1
        self._analyse_case(read_length, prg_structure, expected)

    def test_multipleVariantSites_correctReads(self):
        prg_structure = [
            ['ACA'],
            ['CT', 'GG'],
            ['ACA'],
            ['GC', 'AT'],
        ]
        expected = {'ACT', 'AGG', 'AGC', 'AAT'}
        read_length = 3
        self._analyse_case(read_length, prg_structure, expected)

    def test_readLengthGreaterThanGenome_noReadsGenerated(self):
        prg_structure = [
            ['T'],
            ['CT', 'AG'],
        ]
        expected = set()
        read_length = 4
        self._analyse_case(read_length, prg_structure, expected)

    def test_readLengthGreaterThanGenome_someReadsGenerated(self):
        prg_structure = [
            ['T'],
            ['CT', 'AG'],
            ['A'],
            ['CA', 'AT'],
        ]
        expected = {'TACA', 'TAAT', 'GACA', 'GAAT'}
        read_length = 4
        self._analyse_case(read_length, prg_structure, expected)

    def test_readsOverlappingMultipleVariantSites_correctReads(self):
        prg_structure = [
            ['TA'],
            ['CT', 'GG'],
            ['A'],
            ['GC', 'AT'],
        ]
        expected = {
            'TACT', 'TAGG', 'TAGC',
            'TAAT', 'GAGC', 'GAAT'
        }
        read_length = 4
        self._analyse_case(read_length, prg_structure, expected)

if __name__ == "__main__":
    unittest.main()
