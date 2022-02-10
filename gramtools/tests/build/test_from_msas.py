from unittest import TestCase
from unittest.mock import patch, MagicMock
from typing import List

from pybedtools import BedTool, Interval
from make_prg.from_msa.prg_builder import PrgBuilder
from make_prg.prg_encoder import (
    ENDIANNESS as mp_ENDIANNESS,
    BYTES_PER_INT as mp_BYTES_PER_INT,
)

from gramtools import ENDIANNESS as gram_ENDIANNESS, BYTES_PER_INT as gram_BYTES_PER_INT
from gramtools.commands.common import ints_to_bytes
from gramtools.commands.build.from_msas import (
    PRGAggregator,
    PRGAggregationError,
    PRGDecodeError,
    get_aggregated_prgs,
)


class TestEncodeBinary(TestCase):
    """
    make_prg and vcf_to_prg_string both have code to convert a list of
    integers into a binary file of x bytes per integer with specified endianness.
    (Unfortunately this violates DRY principle)
    Here we make sure they are the same
    """

    def test_binary_encoding_consistency(self):
        self.assertEqual(gram_ENDIANNESS, mp_ENDIANNESS)
        self.assertEqual(gram_BYTES_PER_INT, mp_BYTES_PER_INT)
        self.assertEqual(gram_BYTES_PER_INT, mp_BYTES_PER_INT)


class TestPrgAggregatorIncorrectUses(TestCase):
    def test_translate_non_variant_marker_fails(self):
        # Variant markers are integers >=5
        agg = PRGAggregator()
        with self.assertRaises(PRGAggregationError):
            agg.translate("ref", 4)

    def test_translate_site_marker_more_than_twice_fails(self):
        # Single odd marker is used to encode a variant site
        agg = PRGAggregator()
        agg.translate("ref", 5)
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

    def test_translate_site_marker_twice(self):
        """
        Legacy format support: converts 5A6C5 to 5A6C6
        """
        result = list()
        agg = PRGAggregator()
        result.append(agg.translate("ref", 5))
        result.append(agg.translate("ref", 5))
        expected = [5, 6]
        self.assertEqual(expected, result)

    def test_translate_markers_across_multiple_references(self):
        result = list()
        agg = PRGAggregator()
        for ref, markers in zip(["ref1", "ref2"], [[5, 6, 6]] * 2):
            for marker in markers:
                result.append(agg.translate(ref, marker))
        expected = [5, 6, 6, 7, 8, 8]
        self.assertEqual(expected, result)


def configure_file_mock(file_mock, list_of_prg_ints):
    all_encodings = map(ints_to_bytes, list_of_prg_ints)
    file_mock.return_value = MagicMock()
    file_mock.return_value.read = MagicMock(side_effect=all_encodings)


@patch("builtins.open")
class TestAggregateMultiplePrgs(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.intervals = [
            Interval("ref1", 0, 5, "file1"),
            Interval("ref2", 10, 12, "file2"),
        ]

    def setUp(self):
        self.agg = PRGAggregator()

    def test_invalid_bytes_fails(self, mock_open):
        list_of_prg_ints = [[6, 6]]
        configure_file_mock(mock_open, list_of_prg_ints)
        with self.assertRaises(PRGAggregationError):
            result = get_aggregated_prgs(self.agg, [self.intervals[0]])

    def test_invalid_bytes_fails2(self, mock_open):
        list_of_prg_ints = [[0]]
        configure_file_mock(mock_open, list_of_prg_ints)
        with self.assertRaises(PRGDecodeError):
            result = get_aggregated_prgs(self.agg, [self.intervals[0]])

    def test_single_prg(self, mock_open):
        list_of_prg_ints = [[5, 1, 1, 6, 1, 2, 6]]
        configure_file_mock(mock_open, list_of_prg_ints)
        result = get_aggregated_prgs(self.agg, [self.intervals[0]])
        self.assertEqual(list_of_prg_ints[0], result)

    def test_multiple_prgs(self, mock_open):
        list_of_prg_ints = [[5, 6, 6], [5, 6, 6, 6, 6, 7, 8, 8]]
        configure_file_mock(mock_open, list_of_prg_ints)
        result = get_aggregated_prgs(self.agg, self.intervals)
        expected = [5, 6, 6, 7, 8, 8, 8, 8, 9, 10, 10]
        self.assertEqual(expected, result)
