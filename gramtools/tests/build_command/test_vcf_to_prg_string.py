"""
Tested behaviour:
    * Overlapping vcf records (same CHROM, same POS): the first is kept, the rest is ignored
    * Fasta records with no variation: they are output to the prg string in the order seen in the input fasta.

    * Adjacent records: site markers placed adjacent to each other,
    * Allele markers: in 'normal' mode, use even (allele) marker, in legacy use odd (site) marker at end of site.
"""
import filecmp
from tempfile import mkdtemp
from shutil import rmtree
from pathlib import Path
import os
import subprocess
from unittest import TestCase, mock

from gramtools.commands.build.vcf_to_prg_string import Vcf_to_prg, int_to_bytes
from gramtools.tests.mocks import _MockVcfRecord


"""
In-memory tests, via mocking pysam and fasta loading
"""


@mock.patch("gramtools.commands.build.vcf_to_prg_string.load_fasta")
@mock.patch("gramtools.commands.build.vcf_to_prg_string.VariantFile", spec=True)
class Test_MemberFunctions(TestCase):
    def test_nonACGT_fails_to_convert(self, mock_var_file, mock_load_fasta):
        mock_var_file.return_value.fetch.return_value = iter([])
        mock_load_fasta.return_value = {}
        converter = Vcf_to_prg("", "", "")
        with self.assertRaises(ValueError):
            int_to_bytes("N")


@mock.patch("gramtools.commands.build.vcf_to_prg_string.load_fasta")
@mock.patch("gramtools.commands.build.vcf_to_prg_string.VariantFile", spec=True)
class Test_VcfToPrgString(TestCase):
    chroms = {"ref1": "AGCAGC", "ref2": "CCC", "ref3": "GGG"}

    def test_no_variants_returns_ref_chroms(self, mock_var_file, mock_load_fasta):
        recs = []
        mock_var_file.return_value.fetch.return_value = iter(recs)
        mock_load_fasta.return_value = self.chroms
        converter = Vcf_to_prg("", "", "")
        self.assertEqual("AGCAGCCCCGGG", converter._get_string())

    def test_one_variant_chroms_with_no_vars_in_same_order(
        self, mock_var_file, mock_load_fasta
    ):
        recs = [_MockVcfRecord(2, "G", ["CAAA", "CA"], chrom="ref3")]
        mock_var_file.return_value.fetch.return_value = iter(recs)
        mock_load_fasta.return_value = self.chroms
        converter = Vcf_to_prg("", "", "")
        self.assertEqual("AGCAGCCCCG5G6CAAA6CA6G", converter._get_string())

    def test_two_snps_same_chrom(self, mock_var_file, mock_load_fasta):
        recs = [
            _MockVcfRecord(1, "A", "G", chrom="ref1"),
            _MockVcfRecord(3, "C", ["T", "G"], chrom="ref1"),
        ]
        mock_var_file.return_value.fetch.return_value = iter(recs)
        mock_load_fasta.return_value = self.chroms
        converter = Vcf_to_prg("", "", "")
        self.assertEqual("5A6G6G7C8T8G8AGCCCCGGG", converter._get_string())

    def test_one_ins_and_one_del_diff_chroms(self, mock_var_file, mock_load_fasta):
        recs = [
            _MockVcfRecord(3, "C", ["CGG"], chrom="ref1"),
            _MockVcfRecord(1, "CCC", ["C"], chrom="ref2"),
        ]
        mock_var_file.return_value.fetch.return_value = iter(recs)
        mock_load_fasta.return_value = self.chroms
        converter = Vcf_to_prg("", "", "")
        self.assertEqual("AG5C6CGG6AGC7CCC8C8GGG", converter._get_string())

    def test_adjacent_snps_kept(self, mock_var_file, mock_load_fasta):
        recs = [
            _MockVcfRecord(1, "C", ["G"], chrom="ref2"),
            _MockVcfRecord(2, "C", ["A"], chrom="ref2"),
        ]
        mock_var_file.return_value.fetch.return_value = iter(recs)
        mock_load_fasta.return_value = self.chroms
        converter = Vcf_to_prg("", "", "")
        self.assertEqual("AGCAGC5C6G67C8A8CGGG", converter._get_string())


@mock.patch("gramtools.commands.build.vcf_to_prg_string.load_fasta")
@mock.patch("gramtools.commands.build.vcf_to_prg_string.VariantFile", spec=True)
class Test_Representation(TestCase):
    chroms = {"ref1": "ACACAA"}
    recs = [
        _MockVcfRecord(1, "A", ["G"], chrom="ref1"),
        _MockVcfRecord(5, "A", ["AAA"], chrom="ref1"),
    ]

    def test_legacy_representation(self, mock_var_file, mock_load_fasta):
        mock_var_file.return_value.fetch.return_value = iter(self.recs)
        mock_load_fasta.return_value = self.chroms
        converter = Vcf_to_prg("", "", "", mode="legacy")
        self.assertEqual("5A6G5CAC7A8AAA7A", converter._get_string())

    def test_integer_representation(self, mock_var_file, mock_load_fasta):
        mock_var_file.return_value.fetch.return_value = iter(self.recs)
        mock_load_fasta.return_value = self.chroms
        converter = Vcf_to_prg("", "", "")
        self.assertEqual(
            [5, 1, 6, 3, 6, 2, 1, 2, 7, 1, 8, 1, 1, 1, 8, 1], converter._get_ints()
        )


@mock.patch("gramtools.commands.build.vcf_to_prg_string.load_fasta")
@mock.patch("gramtools.commands.build.vcf_to_prg_string.VariantFile", spec=True)
class Test_VcfToPrgString_OverlappingRecords(TestCase):
    """
    When records overlap, keeps only the first one
    """

    chroms = {"ref1": "TTTT", "ref2": "CCC"}

    def test_snps_at_same_position(self, mock_var_file, mock_load_fasta):
        recs = [
            _MockVcfRecord(1, "TTTT", ["T"], chrom="ref1"),
            _MockVcfRecord(2, "T", ["C"], chrom="ref1"),
        ]
        mock_var_file.return_value.fetch.return_value = iter(recs)
        mock_load_fasta.return_value = self.chroms
        converter = Vcf_to_prg("", "", "")
        self.assertEqual("5TTTT6T6CCC", converter._get_string())

    def test_snp_inside_del(self, mock_var_file, mock_load_fasta):
        recs = [
            _MockVcfRecord(2, "T", ["G"], chrom="ref1"),
            _MockVcfRecord(2, "T", ["C"], chrom="ref1"),
        ]
        mock_var_file.return_value.fetch.return_value = iter(recs)
        mock_load_fasta.return_value = self.chroms
        converter = Vcf_to_prg("", "", "")
        self.assertEqual("T5T6G6TTCCC", converter._get_string())


@mock.patch("gramtools.commands.build.vcf_to_prg_string.load_fasta")
@mock.patch("gramtools.commands.build.vcf_to_prg_string.VariantFile", spec=True)
class Test_Filter_Records(TestCase):
    chroms = {"JAC": "AAC"}

    def test_filter_pass_record_kept(self, mock_var_file, mock_load_fasta):
        recs = [_MockVcfRecord(2, "A", "G", chrom="JAC")]  # Default filter: PASS
        mock_var_file.return_value.fetch.return_value = iter(recs)
        mock_load_fasta.return_value = self.chroms
        converter = Vcf_to_prg("", "", "")
        self.assertEqual(converter._get_string(), "A5A6G6C")

    def test_filter_nonpass_record_skipped(self, mock_var_file, mock_load_fasta):
        recs = [_MockVcfRecord(2, "A", "G", filters={"LOW_QUAL": ""})]
        mock_var_file.return_value.fetch.return_value = iter(recs)
        mock_load_fasta.return_value = self.chroms
        converter = Vcf_to_prg("", "", "")


"""
File-based tests
"""
from gramtools.tests import data_dir

data_dir = data_dir / "vcf_to_prg_string"


class Utility_Tester(object):
    """
    Runs python vcf to prg utility
    """

    def __init__(self, file_prefix):
        self.vcf_file = data_dir / Path(f"{file_prefix}.vcf")
        self.ref_file = data_dir / Path(f"{file_prefix}.ref.fa")
        self.expected_prg = data_dir / Path(f"{file_prefix}.prg")

        self.tmp_dir = mkdtemp()
        # Clumsily below is also the name of the binary output file
        self.outfile_prefix = Path(f"{self.tmp_dir}") / Path("prg")
        self.outfile = Path(str(self.outfile_prefix) + ".prg")

    def _run(self, mode="normal", check=False):
        converter = Vcf_to_prg(
            str(self.vcf_file), str(self.ref_file), self.outfile_prefix, mode=mode
        )
        converter._write_string()

    def _cleanup(self):
        rmtree(self.tmp_dir)


class Test_VcfFormatInputs(TestCase):
    test_cases = {"bad_header_complex_variant"}

    def test_no_fileformat_tag(self):
        """
        Bad header causes pysam error throw
        """
        fname = "bad_header_complex_variant"
        tester = Utility_Tester(fname)

        with self.assertRaises(Exception):
            tester._run(check=True)

    def test_adjacent_snps_allowed(self):
        fname = "two_adjacent_snps_two_alts"
        tester = Utility_Tester(fname)
        tester._run()
        self.assertTrue(filecmp.cmp(tester.expected_prg, tester.outfile, shallow=False))
        tester._cleanup()
