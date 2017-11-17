import filecmp
import os
import subprocess
import unittest

this_file = os.path.abspath(__file__)
vcf_to_linear_prg_pl = os.path.abspath(os.path.join(this_file, os.pardir, os.pardir, 'utils', 'vcf_to_linear_prg.pl'))
data_dir = os.path.abspath(os.path.join(this_file, os.pardir, 'data', 'vcf_to_linear_prg_perl'))

class TestVcfToLinearPrgPerl(unittest.TestCase):
    def test_found_perl_script(self):
        '''Test that we found the script vcf_to_linear_prg.pl'''
        self.assertTrue(os.path.exists(vcf_to_linear_prg_pl))

    # Every test follows the same pattern: run the script and
    # check the output file is the same as the expected output
    # file. The only thing changing here is the names of the files
    def _test_one_run(self, file_prefix):
        vcf_file = os.path.join(data_dir, file_prefix + '.vcf')
        ref_file = os.path.join(data_dir, file_prefix + '.ref.fa')
        expected_prg = os.path.join(data_dir, file_prefix + '.prg')
        outfile = 'tmp.' + file_prefix + '.prg'
        if os.path.exists(outfile):
            os.unlink(outfile)
        command = ' '.join([
            'perl', vcf_to_linear_prg_pl,
            '--vcf', vcf_file,
            '--ref',  ref_file,
            '--min_freq 0.01',
            '--outfile', outfile,
        ])
        completed_process = subprocess.run(command, shell=True)
        self.assertEqual(0, completed_process.returncode)
        self.assertTrue(filecmp.cmp(expected_prg, outfile, shallow=False))
        for suffix in ['', '.fa', '.mask_alleles', '.mask_sites', '.vcf']:
            os.unlink(outfile + suffix)


    def test_prg_with_no_variants(self):
        '''Test make PRG when no variants in VCF'''
        self._test_one_run('prg_with_no_variants')


    def test_prg_with_one_snp(self):
        '''Test make PRG when one SNP VCF'''
        self._test_one_run('prg_with_one_snp')


    def test_prg_with_two_separate_snps(self):
        '''Test make PRG with two separate SNPs'''
        self._test_one_run('prg_with_two_separate_snps')


    def test_prg_with_one_snp_two_alts(self):
        '''Test make PRG with one snp that has two alts'''
        self._test_one_run('prg_with_one_snp_two_alts')


    def test_prg_with_two_adjacent_snps(self):
        '''Test make PRG with two adjacent snps'''
        self._test_one_run('prg_with_two_adjacent_snps')


    def test_prg_with_two_adjacent_snps_two_alts(self):
        '''Test make PRG, two adjacent SNPs, with >1 ALT'''
        self._test_one_run('prg_with_two_adjacent_snps_two_alts')


    def test_prg_with_one_ins_and_one_del(self):
        '''Test prg witn an insertion and deletion independent of each other'''
        self._test_one_run('prg_with_one_ins_and_one_del')


