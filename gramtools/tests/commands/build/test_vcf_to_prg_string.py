"""
Behaviour of conversion utility:
    * Overlapping vcf records (same CHROM, same POS): the first is kept, the rest is ignored
    * Fasta records with no variation: they are output to the prg string,
    concatenated at the end of the prg string, in the order seen in the input fasta.

    * Adjacent records: site markers placed adjacent to each other,
    * Allele markers: in 'normal' mode, the python utility marks the end of a variant site
    with an allele (even) marker, not a site marker.

Bugs/Notes:
        Possible lack of flexibility in vcf input 'bad formatting'
        due to use of PyVcf module.

There used to be a perl utility for conversion. Cases where it worked the same as current utility are
under `test_routine_common_behaviour`
"""
import filecmp
from tempfile import mkdtemp
from shutil import rmtree
from pathlib import Path
import os
import subprocess
import unittest

base_dir = Path(__file__).parent.parent.parent.parent
python_utility = base_dir / "commands" / "build" / "vcf_to_prg_string.py"
data_dir = base_dir / "tests" / "data" / "vcf_to_prg_string"


class Utility_Tester(object):
    """
    Initialised with a file prefix, and able to run perl and python utilities.
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
        """
        :param mode: either "normal" or "legacy". In legacy, variant sites encoded
        in same way as perl utility.
        :return:
        """
        python_command = [
            "python3",
            str(python_utility),
            str(self.vcf_file),
            str(self.ref_file),
            "--outfile",
            self.outfile_prefix,
            "--mode",
            mode,
        ]

        completed_process = subprocess.run(
            python_command, check=check, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        return completed_process

    def _cleanup(self):
        rmtree(self.tmp_dir)


class Test_UtilityScriptsFunction(unittest.TestCase):
    def test_find_utility_script(self):
        """Test that we find the conversion utilities"""
        self.assertTrue(python_utility.exists())


class Test_VcfToPrgString(unittest.TestCase):
    def test_routine_common_behaviour(self):
        test_cases = {
            "no_variants": "Vcf file has no records",
            "one_snp": "",
            "two_separate_snps": "The SNPs are not consecutive",
            "one_snp_two_alts": "Deal with more than one alternatives",
            "one_ins_and_one_del": "Indels",
            "complex_variant": "Mixture of indel & SNP",
            "two_snps_same_place": "Same position covered 1+ times;"
            "Expect only first vcf record kept",
            "snp_inside_del": "Overlapping records;"
            "Expect only first vcf record kept",
        }

        for fname in test_cases.keys():
            with self.subTest(fname=fname):
                tester = Utility_Tester(fname)

                completed_process = tester._run(mode="legacy")
                self.assertEqual(0, completed_process.returncode)
                self.assertTrue(
                    filecmp.cmp(tester.expected_prg, tester.outfile, shallow=False)
                )
                tester._cleanup()

    def test_adjacent_snps_allowed(self):
        test_cases = {
            "two_adjacent_snps": "Two SNPs at consecutive positions.",
            "two_adjacent_snps_two_alts": "Same as above,"
            " with more than one alternative",
        }

        for fname in test_cases.keys():
            with self.subTest(fname=fname):
                tester = Utility_Tester(fname)
                completed_process = tester._run(mode="legacy")
                self.assertEqual(0, completed_process.returncode)
                self.assertTrue(
                    filecmp.cmp(
                        os.path.join(data_dir, "python_" + fname + ".prg"),
                        tester.outfile,
                        shallow=False,
                    )
                )
                tester._cleanup()

    def test_AlleleRepresentation(self):
        """
        Tests the new representation of variant sites implemented
        by default in the python utility (mode="normal")
        """
        test_cases = {
            "one_ins_and_one_del": "",
            "two_adjacent_snps": "Same as above," " with more than one alternative",
        }

        for fname in test_cases.keys():
            with self.subTest(fname=fname):
                tester = Utility_Tester(fname)

                completed_process = tester._run()  # normal mode is default
                self.assertEqual(0, completed_process.returncode)
                self.assertTrue(
                    filecmp.cmp(
                        os.path.join(data_dir, "python_normalMode_" + fname + ".prg"),
                        tester.outfile,
                        shallow=False,
                    )
                )
                tester._cleanup()


class Test_VcfFormatInputs(unittest.TestCase):
    test_cases = {"bad_header_complex_variant"}

    def test_singlehash_headers(self):
        """
        Note/BUG here: it is probably a PyVcf limitation that it does not parse
        this file correctly.
        """
        fname = "bad_header_complex_variant"
        tester = Utility_Tester(fname)

        with self.assertRaises(Exception):
            tester._run(check=True)
