import os
import unittest

this_file = os.path.abspath(__file__)
vcf_to_linear_prg_pl = os.path.abspath(os.path.join(this_file, os.pardir, os.pardir, 'utils', 'vcf_to_linear_prg.pl'))

class TestVcfToLinearPrgPerl(unittest.TestCase):
    def test_found_perl_script(self):
        '''Test that we found the script vcf_to_linear_prg.pl'''
        self.assertTrue(os.path.exists(vcf_to_linear_prg_pl))


