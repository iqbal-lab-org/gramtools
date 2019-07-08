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
data_dir = os.path.abspath(os.path.join(this_file, os.pardir, 'data', 'vcf_to_linear_prg_perl'))


class Test_VcfToPrgString(unittest.TestCase):
    def test_found_perl_script(self):
        """Test that we found the conversion utilities"""
        self.assertTrue(os.path.exists(perl_utility))
        self.assertTrue(os.path.exists(python_utility))

    # Every test follows the same pattern: run the script and
    # check the output file is the same as the expected output
    # file. The only thing changing here is the names of the files
    def _test_one_run(self, file_prefix):
        vcf_file = os.path.join(data_dir, file_prefix + '.vcf')
        ref_file = os.path.join(data_dir, file_prefix + '.ref.fa')
        expected_prg = os.path.join(data_dir, file_prefix + '.prg')
        expected_vcf = os.path.join(data_dir, file_prefix + '.expect.vcf')
        outfile = 'tmp.' + file_prefix + '.prg'
        if os.path.exists(outfile):
            os.unlink(outfile)

        # Test the perl utility
        perl_command = ' '.join([
            'perl', perl_utility,
            '--vcf', vcf_file,
            '--ref', ref_file,
            '--min_freq 0.01',
            '--outfile', outfile,
        ])
        completed_process = subprocess.run(perl_command, shell=True)
        self.assertEqual(0, completed_process.returncode)
        self.assertTrue(filecmp.cmp(expected_prg, outfile, shallow=False))
        self.assertTrue(filecmp.cmp(expected_vcf, outfile + '.vcf', shallow=False))
        for suffix in ['', '.fa', '.mask_alleles', '.mask_sites', '.vcf']:
            os.unlink(outfile + suffix)

        python_command = ' '.join([
            'python3', python_utility,
            vcf_file,
            ref_file,
            '--outfile', outfile,
            '--mode', "legacy",
        ])
        completed_process = subprocess.run(python_command, shell=True)
        self.assertEqual(0, completed_process.returncode)
        self.assertTrue(filecmp.cmp(expected_prg, outfile, shallow=False))

        os.unlink(outfile)


    def test_prg_with_no_variants(self):
        """Test make PRG when no variants in VCF"""
        self._test_one_run('prg_with_no_variants')

    def test_prg_with_one_snp(self):
        """Test make PRG when one SNP VCF"""
        self._test_one_run('prg_with_one_snp')

    def test_prg_with_two_separate_snps(self):
        """Test make PRG with two separate SNPs"""
        self._test_one_run('prg_with_two_separate_snps')

    def test_prg_with_one_snp_two_alts(self):
        """Test make PRG with one snp that has two alts"""
        self._test_one_run('prg_with_one_snp_two_alts')

    def test_prg_with_two_adjacent_snps(self):
        """Test make PRG with two adjacent snps"""
        self._test_one_run('prg_with_two_adjacent_snps')

    def test_prg_with_two_adjacent_snps_two_alts(self):
        """Test make PRG, two adjacent SNPs, with >1 ALT"""
        self._test_one_run('prg_with_two_adjacent_snps_two_alts')

    def test_prg_with_one_ins_and_one_del(self):
        """Test prg witn an insertion and deletion independent of each other"""
        self._test_one_run('prg_with_one_ins_and_one_del')

    def test_prg_with_complex_variant(self):
        """Test PRG when we have complex variant"""
        self._test_one_run('prg_with_complex_variant')

    # The current behaviour of two SNPs in the same place but in separate
    # VCF lines is to only use the first one.  Testing this here to
    # document the behaviour as much as anything else.
    # We may want to fix this in the future.
    def test_prg_with_two_snps_same_place(self):
        """Test PRG with two SNPs in the same plave in separate VCF lines"""
        self._test_one_run('prg_with_two_snps_same_place')

    # If there is one variant inside another one, eg this test is:
    # Ref:  AGCACGTCATGA
    # var1: AG-----CATGA
    # var2: AGCTCGTCATGA
    #          * <- SNP here!
    # Then only var1 is put into the graph, but both lines
    # of VCF are put into the output.
    # This should be fixed in the future
    def test_prg_with_snp_inside_del(self):
        """Test when snp inside deletion"""
        self._test_one_run('prg_with_snp_inside_del')
