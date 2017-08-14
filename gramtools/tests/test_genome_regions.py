import unittest

from . import common
from .. import prg


class TestDirectionalRegionRange(unittest.TestCase):
    def test_givenPrg_regionAllelesSetCorrectly(self):
        prg_structure = [
            ['C'],
            ['G', 'T'],
            ['A'],
            ['TC', 'A'],
        ]
        expected = [
            ['C'],
            ['G', 'T'],
            ['A'],
            ['TC', 'A'],
        ]
        prg_seq = common.compose_prg(prg_structure)
        regions = prg.parse(prg_seq)
        result = []
        for region in regions:
            result.append([''.join(x) for x in region.alleles])
        self.assertEqual(result, expected)

    def test_reverseRange_yieldRegionsInReverse(self):
        prg_structure = [
            ['C'],
            ['G', 'T'],
            ['A'],
            ['TC', 'A'],
        ]
        expected = [
            ['A'],
            ['G', 'T'],
            ['C'],
        ]
        prg_seq = common.compose_prg(prg_structure)
        regions = prg.parse(prg_seq)

        start_region = regions[-1]
        regions_range = regions.range(start_region, reverse=True)

        result = []
        for region in regions_range:
            result.append([''.join(x) for x in region.alleles])
        self.assertEqual(result, expected)

    def test_forwardRange_yieldRegionsInOrder(self):
        prg_structure = [
            ['C'],
            ['G', 'T'],
            ['A'],
            ['TC', 'A'],
        ]
        expected = [
            ['G', 'T'],
            ['A'],
            ['TC', 'A'],
        ]
        prg_seq = common.compose_prg(prg_structure)
        regions = prg.parse(prg_seq)

        start_region = regions[0]
        directional_regions = list(regions.range(start_region, reverse=False))
        results = []
        for region in directional_regions:
            results.append([''.join(x) for x in region.alleles])

        self.assertEqual(results, expected)


    def test_reverseRange_startRegionNotInRange(self):
        prg_structure = [
            ['C'],
            ['G', 'T'],
            ['A'],
            ['TC', 'A'],
        ]
        prg_seq = common.compose_prg(prg_structure)
        regions = prg.parse(prg_seq)

        start_region = regions[-1]
        range_alleles = [r.alleles for r in
                         regions.range(start_region, reverse=True)]
        self.assertNotIn(start_region.alleles, range_alleles)

    def test_forwardRange_startRegionNotInRange(self):
        prg_structure = [
            ['C'],
            ['G', 'T'],
            ['A'],
            ['TC', 'A'],
        ]
        prg_seq = common.compose_prg(prg_structure)
        regions = prg.parse(prg_seq)

        start_region = regions[1]
        range_alleles = [r.alleles for r in
                         regions.range(start_region, reverse=False)]
        self.assertNotIn(start_region.alleles, range_alleles)
