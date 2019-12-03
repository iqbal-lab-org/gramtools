import unittest
from gramtools.commands.simulate import prg_regions_parser


class TestParsePrg(unittest.TestCase):
    def test_noVariantSites_oneBlock(self):
        prg_seq = "ACTGTCCTC"
        regions = prg_regions_parser.parse(prg_seq)
        regions = [region for region in regions]
        self.assertEqual(len(regions), 1)

    def test_noVariantSites_oneAlleleBlock(self):
        prg_seq = "ACTGTCCTC"
        regions = prg_regions_parser.parse(prg_seq)
        region = [region for region in regions][0]
        result = ["".join(x) for x in region.alleles]
        expected = ["ACTGTCCTC"]
        self.assertEqual(result, expected)

    def test_twoVariantSites_correctNumVariantSites(self):
        prg_seq = "AC5A6T6A5CC7CA8T7CC"
        regions = prg_regions_parser.parse(prg_seq)

        num_variant_sites = sum(region.is_variant_site for region in regions)
        self.assertEqual(num_variant_sites, 2)

    def test_twoVariantSites_correctVariantMarkers(self):
        prg_seq = "ACTGT5A6T6A5CC7CA8T7CC"
        regions = prg_regions_parser.parse(prg_seq)

        variant_markers = []
        for region in regions:
            if region.is_variant_site:
                variant_markers.append(region.variant_site_marker)

        expected = ["5", "7"]
        self.assertEqual(variant_markers, expected)

    def test_nonVariantSites_singleAllelePerBlock(self):
        prg_seq = "ACTGT5A6T6A5CC7CA8T7CC"
        regions = prg_regions_parser.parse(prg_seq)

        for region in regions:
            if not region.is_variant_site:
                self.assertEqual(len(region.alleles), 1)

    def test_threeVariantAlleles_orderedInAlleleAttr(self):
        prg_seq = "TT5A6T6AA5CC"
        regions = prg_regions_parser.parse(prg_seq)
        variant_site_block = [region for region in regions if region.is_variant_site][0]
        expected = ["A", "T", "AA"]
        result = ["".join(x) for x in variant_site_block.alleles]
        self.assertEqual(result, expected)

    def test_rightEdgeVariantSite_correctBlocks(self):
        prg_seq = "TT5A6T5AA7C8A7"
        regions = prg_regions_parser.parse(prg_seq)

        expected = [
            ["TT"],
            ["A", "T"],  # first variant site
            ["AA"],
            ["C", "A"],  # second variant site
        ]

        for region, expected in zip(regions, expected):
            result = ["".join(x) for x in region.alleles]
            self.assertEqual(result, expected)

    def test_leftEdgeVariantSite_correctBlocks(self):
        prg_seq = "5A6T5AA7C8A7CC"
        regions = prg_regions_parser.parse(prg_seq)

        expected = [
            ["A", "T"],  # first variant site
            ["AA"],
            ["C", "A"],  # second variant site
            ["CC"],
        ]

        for region, expected in zip(regions, expected):
            result = ["".join(x) for x in region.alleles]
            self.assertEqual(result, expected)

    def test_bothEdgesVariantSite_correctBlocks(self):
        prg_seq = "5A6T5AA7C8A7"
        regions = prg_regions_parser.parse(prg_seq)

        expected = [
            ["A", "T"],  # first variant site
            ["AA"],
            ["C", "A"],  # second variant site
        ]

        for region, expected in zip(regions, expected):
            result = ["".join(x) for x in region.alleles]
            self.assertEqual(result, expected)

    def test_singleVariantSite_singleVariantBlock(self):
        prg_seq = "5A6T5"
        regions = prg_regions_parser.parse(prg_seq)

        expected = [["A", "T"]]

        for region, expected in zip(regions, expected):
            result = ["".join(x) for x in region.alleles]
            self.assertEqual(result, expected)

    def test_multiCharVariantSiteMarker_correctBlocks(self):
        prg_seq = "TT97AC98GT97CC"
        regions = prg_regions_parser.parse(prg_seq)

        expected = [["TT"], ["AC", "GT"], ["CC"]]

        for region, expected in zip(regions, expected):
            result = ["".join(x) for x in region.alleles]
            self.assertEqual(result, expected)

    def test_leftEdgeMultiCharVarSiteMarker_correctBlocks(self):
        prg_seq = "97AC98GT97CC"
        regions = prg_regions_parser.parse(prg_seq)

        expected = [["AC", "GT"], ["CC"]]

        for region, expected in zip(regions, expected):
            result = ["".join(x) for x in region.alleles]
            self.assertEqual(result, expected)

    def test_rightEdgeMultiCharVarSiteMarker_correctBlocks(self):
        prg_seq = "TT97AC98GT97"
        regions = prg_regions_parser.parse(prg_seq)

        expected = [["TT"], ["AC", "GT"]]

        for region, expected in zip(regions, expected):
            result = ["".join(x) for x in region.alleles]
            self.assertEqual(result, expected)

    def test_bothEdgesMultiCharVarSiteMarker_correctBlocks(self):
        prg_seq = "97AC98GT97"
        regions = prg_regions_parser.parse(prg_seq)

        expected = [["AC", "GT"]]

        for region, expected in zip(regions, expected):
            result = ["".join(x) for x in region.alleles]
            self.assertEqual(result, expected)


class TestIterPeek(unittest.TestCase):
    def test_evenNumElements_correctPeekValues(self):
        prg_seq = "ACTG"
        iter_prg = prg_regions_parser.IterPeek(prg_seq)

        expected = ["C", "T", "G", None]
        peek_values = [iter_prg.peek for _ in iter_prg]
        self.assertEqual(peek_values, expected)

    def test_oddNumElements_correctPeekValues(self):
        prg_seq = "ACTGA"
        iter_prg = prg_regions_parser.IterPeek(prg_seq)

        expected = ["C", "T", "G", "A", None]
        peek_values = [iter_prg.peek for _ in iter_prg]
        self.assertEqual(peek_values, expected)

    def test_singleElement_nonePeek(self):
        prg_seq = "A"
        iter_prg = prg_regions_parser.IterPeek(prg_seq)

        expected = [None]
        peek_values = [iter_prg.peek for _ in iter_prg]
        self.assertEqual(peek_values, expected)

    def test_singleElement_correctForloop(self):
        prg_seq = "A"
        iter_prg = prg_regions_parser.IterPeek(prg_seq)

        expected = ["A"]
        values = [x for x in iter_prg]
        self.assertEqual(values, expected)

    def test_evenElement_correctForloop(self):
        prg_seq = "ACTG"
        iter_prg = prg_regions_parser.IterPeek(prg_seq)

        expected = ["A", "C", "T", "G"]
        values = [x for x in iter_prg]
        self.assertEqual(values, expected)

    def test_oddElement_correctForloop(self):
        prg_seq = "ACTGA"
        iter_prg = prg_regions_parser.IterPeek(prg_seq)

        expected = ["A", "C", "T", "G", "A"]
        values = [x for x in iter_prg]
        self.assertEqual(values, expected)


if __name__ == "__main__":
    unittest.main()
