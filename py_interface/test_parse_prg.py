import unittest
import parse_prg


class TestParsePrg(unittest.TestCase):
    def test_noVariantSites_oneBlock(self):
        prg = 'ACTGTCCTC'
        regions = parse_prg.parse(prg)
        regions = [region for region in regions]
        self.assertEqual(len(regions), 1)

    def test_noVariantSites_oneAlleleBlock(self):
        prg = 'ACTGTCCTC'
        regions = parse_prg.parse(prg)
        region = [region for region in regions][0]
        self.assertEqual(region.alleles, ['ACTGTCCTC'])

    def test_twoVariantSites_correctNumVariantSites(self):
        prg = 'AC5A6T6A5CC7CA8T7CC'
        regions = parse_prg.parse(prg)

        num_variant_sites = sum(region.is_variant_site
                                for region in regions)
        self.assertEqual(num_variant_sites, 2)

    def test_twoVariantSites_correctVariantMarkers(self):
        prg = 'ACTGT5A6T6A5CC7CA8T7CC'
        regions = parse_prg.parse(prg)

        variant_markers = []
        for region in regions:
            if region.is_variant_site:
                variant_markers.append(region.variant_site_marker)

        expected = ['5', '7']
        self.assertEqual(variant_markers, expected)

    def test_nonVariantSites_singleAllelePerBlock(self):
        prg = 'ACTGT5A6T6A5CC7CA8T7CC'
        regions = parse_prg.parse(prg)

        for region in regions:
            if not region.is_variant_site:
                self.assertEqual(len(region.alleles), 1)

    def test_threeVariantAlleles_orderedInAlleleAttr(self):
        prg = 'TT5A6T6AA5CC'
        regions = parse_prg.parse(prg)
        variant_site_block = [region for region in regions
                              if region.is_variant_site][0]
        expected = ['A', 'T', 'AA']
        self.assertEqual(variant_site_block.alleles, expected)

    def test_rightEdgeVariantSite_correctBlocks(self):
        prg = 'TT5A6T5AA7C8A7'
        regions = parse_prg.parse(prg)

        expected = [
            ['TT'],
            ['A', 'T'],  # first variant site
            ['AA'],
            ['C', 'A'],  # second variant site
        ]

        for region, expected in zip(regions, expected):
            self.assertEqual(region.alleles, expected)

    def test_leftEdgeVariantSite_correctBlocks(self):
        prg = '5A6T5AA7C8A7CC'
        regions = parse_prg.parse(prg)

        expected = [
            ['A', 'T'],  # first variant site
            ['AA'],
            ['C', 'A'],  # second variant site
            ['CC']
        ]

        for region, expected in zip(regions, expected):
            self.assertEqual(region.alleles, expected)

    def test_bothEdgesVariantSite_correctBlocks(self):
        prg = '5A6T5AA7C8A7'
        regions = parse_prg.parse(prg)

        expected = [
            ['A', 'T'],  # first variant site
            ['AA'],
            ['C', 'A'],  # second variant site
        ]

        for region, expected in zip(regions, expected):
            self.assertEqual(region.alleles, expected)

    def test_singleVariantSite_singleVariantBlock(self):
        prg = '5A6T5'
        regions = parse_prg.parse(prg)

        expected = [
            ['A', 'T'],
        ]

        for region, expected in zip(regions, expected):
            self.assertEqual(region.alleles, expected)

    def test_multiCharVariantSiteMarker_correctBlocks(self):
        prg = 'TT97AC98GT97CC'
        regions = parse_prg.parse(prg)

        expected = [
            ['TT'],
            ['AC', 'GT'],
            ['CC'],
        ]

        for region, expected in zip(regions, expected):
            self.assertEqual(region.alleles, expected)

    def test_leftEdgeMultiCharVarSiteMarker_correctBlocks(self):
        prg = '97AC98GT97CC'
        regions = parse_prg.parse(prg)

        expected = [
            ['AC', 'GT'],
            ['CC'],
        ]

        for region, expected in zip(regions, expected):
            self.assertEqual(region.alleles, expected)

    def test_rightEdgeMultiCharVarSiteMarker_correctBlocks(self):
        prg = 'TT97AC98GT97'
        regions = parse_prg.parse(prg)

        expected = [
            ['TT'],
            ['AC', 'GT'],
        ]

        for region, expected in zip(regions, expected):
            self.assertEqual(region.alleles, expected)

    def test_bothEdgesMultiCharVarSiteMarker_correctBlocks(self):
        prg = '97AC98GT97'
        regions = parse_prg.parse(prg)

        expected = [
            ['AC', 'GT'],
        ]

        for region, expected in zip(regions, expected):
            self.assertEqual(region.alleles, expected)


class TestIterPeek(unittest.TestCase):
    def test_evenNumElements_correctPeekValues(self):
        prg = 'ACTG'
        iter_prg = parse_prg.IterPeek(prg)

        expected = ['C', 'T', 'G', None]
        peek_values = [iter_prg.peek for _ in iter_prg]
        self.assertEqual(peek_values, expected)

    def test_oddNumElements_correctPeekValues(self):
        prg = 'ACTGA'
        iter_prg = parse_prg.IterPeek(prg)

        expected = ['C', 'T', 'G', 'A', None]
        peek_values = [iter_prg.peek for _ in iter_prg]
        self.assertEqual(peek_values, expected)

    def test_singleElement_nonePeek(self):
        prg = 'A'
        iter_prg = parse_prg.IterPeek(prg)

        expected = [None]
        peek_values = [iter_prg.peek for _ in iter_prg]
        self.assertEqual(peek_values, expected)

    def test_singleElement_correctForloop(self):
        prg = 'A'
        iter_prg = parse_prg.IterPeek(prg)

        expected = ['A']
        values = [x for x in iter_prg]
        self.assertEqual(values, expected)

    def test_evenElement_correctForloop(self):
        prg = 'ACTG'
        iter_prg = parse_prg.IterPeek(prg)

        expected = ['A', 'C', 'T', 'G']
        values = [x for x in iter_prg]
        self.assertEqual(values, expected)

    def test_oddElement_correctForloop(self):
        prg = 'ACTGA'
        iter_prg = parse_prg.IterPeek(prg)

        expected = ['A', 'C', 'T', 'G', 'A']
        values = [x for x in iter_prg]
        self.assertEqual(values, expected)


if __name__ == '__main__':
    unittest.main()
