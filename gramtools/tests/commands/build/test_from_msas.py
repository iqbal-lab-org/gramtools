from unittest import TestCase

from gramtools.commands.build.from_msas import PRGAggregator, PRGAggregationError


class TestPrgAggregatorIncorrectUses(TestCase):
    def test_translate_non_variant_marker_fails(self):
        # Variant markers are integers >=5
        agg = PRGAggregator()
        with self.assertRaises(PRGAggregationError):
            agg.translate("ref", 4)

    def test_translate_site_marker_multiple_times_fails(self):
        # Single odd marker is used to encode a variant site
        agg = PRGAggregator()
        agg.translate("ref", 5)
        with self.assertRaises(PRGAggregationError):
            agg.translate("ref", 5)

    def test_translate_allele_marker_without_site_marker_fails(self):
        # Allele marker comes after site marker in encoding
        agg = PRGAggregator()
        with self.assertRaises(PRGAggregationError):
            agg.translate("ref", 6)


class TestPrgAggregatorCorrectUses(TestCase):
    def test_first_allocated_marker_is_fixed(self):
        # This should not happen in practice, first encountered
        # marker should be 5, but is not this class's job to enforce
        agg = PRGAggregator()
        self.assertEqual(5, agg.translate("ref", 101))

    def test_translate_site_then_allele_marker(self):
        result = list()
        agg = PRGAggregator()
        expected = [5, 6, 6]
        for marker in expected:
            result.append(agg.translate("ref", marker))
        self.assertEqual(expected, result)

    def test_translate_markers_across_multiple_references(self):
        result = list()
        agg = PRGAggregator()
        for ref, markers in zip(["ref1", "ref2"], [[5, 6, 6]] * 2):
            for marker in markers:
                result.append(agg.translate(ref, marker))
        expected = [5, 6, 6, 7, 8, 8]
        self.assertEqual(expected, result)
