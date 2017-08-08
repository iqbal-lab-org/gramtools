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

        for region, expected_alleles in zip(regions, expected):
            self.assertEqual(region.alleles, expected_alleles)

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

        for region, expected_alleles in zip(regions_range, expected):
            self.assertEqual(region.alleles, expected_alleles)

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
        regions_range = regions.range(start_region, reverse=False)

        for region, expected_alleles in zip(regions_range, expected):
            self.assertEqual(region.alleles, expected_alleles)

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
