"""
Two utilities are tested here for making a PRG string from a vcf.
They are expected to have the same behaviour for:
    * Overlapping vcf records (same CHROM, same POS): the first is kept, the rest is ignored
    * Fasta records with no variation: they are output to the prg string. For python utility
     they are concatenated at the end of the prg string, in the order seen in the input fasta.

They have DIFFERENT behaviour for:
    * File production: the python utility only produces a prg string file
    * Adjacent records: the python utility places site markers adjacent to each other,
    while the perl one enumerates all combinations.
    * Allele markers: in 'normal' mode, the python utility marks the end of a variant site
    with an allele (even) marker, not a site marker.

Bugs/Notes:
    Perl utility:
        If there is one variant inside another one, eg this test is:
        Ref:  AGCACGTCATGA
        var1: AG-----CATGA
        var2: AGCTCGTCATGA
                 * <- SNP here!
        Then only var1 is put into the graph, but both lines
        of VCF are put into the output.

    Python utility:
        Possible lack of flexibility in vcf input 'bad formatting'
        due to use of PyVcf module.


"""
import filecmp
import os
import subprocess
import unittest

this_file = os.path.abspath(__file__)
perl_utility = os.path.abspath(os.path.join(this_file, os.pardir,
                                            os.pardir, 'utils',
                                                    'vcf_to_linear_prg.pl'))
python_utility = os.path.abspath(os.path.join(this_file, os.pardir,
                                            os.pardir, 'utils',
                                            'vcf_to_prg_string.py'))
data_dir = os.path.abspath(os.path.join(this_file, os.pardir, 'data', 'vcf_to_prg_string'))

class Utility_Tester(object):
    """
    Initialised with a file prefix, and able to run perl and python utilities.
    """
    def __init__(self, file_prefix):
        self.vcf_file = os.path.join(data_dir, file_prefix + '.vcf')
        self.ref_file = os.path.join(data_dir, file_prefix + '.ref.fa')
        self.expected_prg = os.path.join(data_dir, file_prefix + '.prg')
        self.expected_vcf = os.path.join(data_dir, file_prefix + '.expect.vcf') # For perl util only
        self.outfile_prefix = 'tmp.' + file_prefix
        self.outfile = self.outfile_prefix + '.prg'

    def _run_perl_utility(self):
        if os.path.exists(self.outfile):
            os.unlink(self.outfile)

        # Test the perl utility
        perl_command = ' '.join([
            'perl', perl_utility,
            '--vcf', self.vcf_file,
            '--ref', self.ref_file,
            '--min_freq 0.01',
            '--outfile', self.outfile,
        ])
        completed_process = subprocess.run(perl_command, shell=True)
        return completed_process

    def _perl_cleanup(self):
        for suffix in ['', '.fa', '.mask_alleles', '.mask_sites', '.vcf']:
            os.unlink(self.outfile + suffix)

    def _run_python_utility(self, mode="normal", check=False):
        """
        :param mode: either "normal" or "legacy". In legacy, variant sites encoded
        in same way as perl utility.
        :return:
        """
        python_command = ' '.join([
            'python3', python_utility,
            self.vcf_file,
            self.ref_file,
            '--outfile', self.outfile_prefix,
            '--mode', mode,
        ])

        completed_process = subprocess.run(python_command, shell=True,
                                           check=check)
        return completed_process

    def _python_cleanup(self):
        os.unlink(self.outfile_prefix)
        os.unlink(f'{self.outfile}')


class Test_UtilityScriptsFunction(unittest.TestCase):
    def test_find_utility_script(self):
        """Test that we find the conversion utilities"""
        self.assertTrue(os.path.exists(perl_utility))
        self.assertTrue(os.path.exists(python_utility))


class Test_VcfToPrgString(unittest.TestCase):

    def test_routine_common_behaviour(self):
        """
        This set runs tests where we expect the two utilities
        to produce the exact same outputs.
        """
        # Here's a description of the test cases. There are files that correspond to each.
        test_cases = {'no_variants': "Vcf file has no records",
                      'one_snp': "",
                      'two_separate_snps': "The SNPs are not consecutive",
                      'one_snp_two_alts': "Deal with more than one alternatives",
                      'one_ins_and_one_del': "Indels",
                      'complex_variant': "Mixture of indel & SNP",
                      'two_snps_same_place': "Same position covered 1+ times;"
                                                      "Expect only first vcf record kept",
                      'snp_inside_del': "Overlapping records;"
                                                 "Expect only first vcf record kept"}

        for fname in test_cases.keys():
            with self.subTest(fname=fname):
                tester = Utility_Tester(fname)
                completed_process = tester._run_perl_utility()
                self.assertEqual(0, completed_process.returncode)
                self.assertTrue(filecmp.cmp(tester.expected_prg, tester.outfile, shallow=False))
                self.assertTrue(filecmp.cmp(tester.expected_vcf, tester.outfile + '.vcf', shallow=False))
                tester._perl_cleanup()

                completed_process = tester._run_python_utility(mode="legacy")
                self.assertEqual(0, completed_process.returncode)
                self.assertTrue(filecmp.cmp(tester.expected_prg, tester.outfile, shallow=False))
                tester._python_cleanup()

    def test_routine_unique_behaviour(self):
        """
        This set runs tests where we expect the two utilities
        to behave differently
        """
        test_cases = {
            'two_adjacent_snps': "Two SNPs at consecutive positions."
                                          "The perl utility enumerates all combinations"
                                          "which avoids adjacent site markers in the prg string."
                                          "The python utility lets adjacent site markers"
                                          "occur.",
            'two_adjacent_snps_two_alts': "Same as above,"
                                                   " with more than one alternative",
        }

        for fname in test_cases.keys():
            with self.subTest(fname=fname):
                tester = Utility_Tester(fname)
                completed_process = tester._run_perl_utility()
                self.assertEqual(0, completed_process.returncode)
                self.assertTrue(filecmp.cmp(os.path.join(data_dir, "perl_"+ fname+".prg"),
                                            tester.outfile, shallow=False))
                self.assertTrue(filecmp.cmp(tester.expected_vcf,
                                            tester.outfile + '.vcf', shallow=False))
                tester._perl_cleanup()

                completed_process = tester._run_python_utility(mode="legacy")
                self.assertEqual(0, completed_process.returncode)
                self.assertTrue(filecmp.cmp(os.path.join(data_dir, "python_"+fname+".prg"),
                                            tester.outfile, shallow=False))
                tester._python_cleanup()

    def test_routine_python_normalMode(self):
        """
        Tests the new representation of variant sites implemented
        by default in the python utility (mode="normal")
        """
        test_cases = {
            'one_ins_and_one_del': "" ,
            'two_adjacent_snps': "Same as above,"
                                          " with more than one alternative",
        }

        for fname in test_cases.keys():
            with self.subTest(fname=fname):
                tester = Utility_Tester(fname)

                completed_process = tester._run_python_utility() # No mode is normal mode
                self.assertEqual(0, completed_process.returncode)
                self.assertTrue(filecmp.cmp(os.path.join(data_dir,
                                                         "python_normalMode_"+fname+".prg"),
                                            tester.outfile, shallow=False))
                tester._python_cleanup()


class Test_VcfFormatInputs(unittest.TestCase):
    """
    A note about assertRaises: if ran using the UtilityTester module's functions,
    assertRaises checks whether python subprocess module raises an exception (not
    the python conversion utility itself).

    This is enforced by asking subprocess.run to check its return code.
    """
    test_cases = {'bad_header_complex_variant'}
    def test_singlehash_headers(self):
        """
        Note/BUG here: it is probably a PyVcf limitation that it does not parse
        this file correctly.
        """
        fname='bad_header_complex_variant'
        tester = Utility_Tester(fname)

        with self.assertRaises(Exception):
            tester._run_python_utility(check=True)
