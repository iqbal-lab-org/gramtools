import os
import unittest

from minos import genotyper

modules_dir = os.path.dirname(os.path.abspath(genotyper.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'genotyper')

class TestGenotyper(unittest.TestCase):
    '''test _singleton_alleles_and_coverage'''
    def test_singleton_alleles_and_coverage(self):
        allele_combination_cov = {'1': 20, '3': 1}
        allele_groups_dict = {'1': {0}, '2': {1}, '3': {1,2}, '4': {5,6}}
        self.assertEqual({0: 20}, genotyper.Genotyper._singleton_alleles_and_coverage(allele_combination_cov, allele_groups_dict))
        allele_combination_cov['2'] = 42
        self.assertEqual({0: 20, 1: 42}, genotyper.Genotyper._singleton_alleles_and_coverage(allele_combination_cov, allele_groups_dict))


    def test_total_coverage(self):
        '''test _total_coverage'''
        self.assertEqual(0, genotyper.Genotyper._total_coverage({}))
        self.assertEqual(1, genotyper.Genotyper._total_coverage({'x': 1}))
        self.assertEqual(42, genotyper.Genotyper._total_coverage({'x': 1, 'y': 41}))


    def test_coverage_of_one_haploid_allele(self):
        '''test _coverage_of_one_haploid_allele'''
        allele_combination_cov = {'1': 20, '2': 1}
        allele_groups_dict = {'0': {0}, '1': {1}, '2': {1,2}, '3': {5,6}}
        self.assertEqual(0, genotyper.Genotyper._coverage_of_one_haploid_allele(0, allele_combination_cov, allele_groups_dict))
        self.assertEqual(21, genotyper.Genotyper._coverage_of_one_haploid_allele(1, allele_combination_cov, allele_groups_dict))


    def test_coverage_of_diploid_alleles(self):
        'test _coverage_of_diploid_alleles'''
        # This is like:
        #   {0}: 20
        #   {1}: 80
        #   {0,1}: 10
        #   {0,2): 3
        allele_combination_cov = {'1': 20, '2': 80, '3': 10, '4': 3}
        allele_groups_dict = {'1': {0}, '2': {1}, '3': {0,1}, '4': {0,2}}
        singleton_alleles_cov = genotyper.Genotyper._singleton_alleles_and_coverage(allele_combination_cov, allele_groups_dict)
        self.assertEqual((25, 88), genotyper.Genotyper._coverage_of_diploid_alleles(0, 1, allele_combination_cov, allele_groups_dict, singleton_alleles_cov))


    def test_log_likelihood_homozygous(self):
        '''test _log_likelihood_homozygous'''
        self.assertEqual(-26.71, round(genotyper.Genotyper._log_likelihood_homozygous(100, 90, 95, 0.01, 5, 5), 2))
        self.assertEqual(-44.54, round(genotyper.Genotyper._log_likelihood_homozygous(10, 1, 9, 0.01, 5, 5), 2))


    def test_log_likelihood_heterozygous(self):
        '''test _log_likelihood_heterozygous'''
        #_log_likelihood_heterozygous(cls, mean_depth, allele_depth1, allele_depth2, total_depth,
        #    error_rate, allele_length1, allele_length2, non_zeros1, non_zeros2):
        self.assertEqual(-52.97, round(genotyper.Genotyper._log_likelihood_heterozygous(100, 45, 40, 95, 0.01, 3, 3, 3, 3), 2))
        self.assertEqual(-152.97, round(genotyper.Genotyper._log_likelihood_heterozygous(100, 45, 40, 95, 0.01, 3, 3, 2, 2), 2))


    def test_calculate_log_likelihoods(self):
        '''test _calculate_log_likelihoods'''
        mean_depth = 20
        error_rate = 0.01
        allele_combination_cov = {'1': 2, '2': 20, '3':1}
        allele_groups_dict = {'1': {0}, '2': {1}, '3': {0,1}, '4': {2}}
        allele_per_base_cov = [[0, 1], [20, 19]]
        gtyper = genotyper.Genotyper(mean_depth, error_rate, allele_combination_cov, allele_per_base_cov, allele_groups_dict)
        gtyper._calculate_log_likelihoods()
        expected = [
            ({1}, -11.68),
            ({0, 1}, -26.98),
            ({0}, -124.91),
        ]
        self.assertEqual(3, len(gtyper.likelihoods))
        gtyper.likelihoods = [(x[0], round(x[1], 2)) for x in gtyper.likelihoods]
        self.assertEqual(expected, gtyper.likelihoods)


    def test_run(self):
        '''test run'''
        mean_depth = 20
        error_rate = 0.01
        allele_combination_cov = {'1': 2, '2': 20, '3':1}
        allele_groups_dict = {'1': {0}, '2': {1}, '3': {0,1}, '4': {2}}
        allele_per_base_cov = [[0, 1], [20, 19]]
        gtyper = genotyper.Genotyper(mean_depth, error_rate, allele_combination_cov, allele_per_base_cov, allele_groups_dict)
        expected = [
            ({1}, -11.68),
            ({0, 1}, -26.98),
            ({0}, -124.91),
        ]
        gtyper.run()
        self.assertEqual(len(expected), len(gtyper.likelihoods))
        for i in range(len(expected)):
            self.assertEqual(expected[i][0], gtyper.likelihoods[i][0])
            self.assertAlmostEqual(expected[i][1], gtyper.likelihoods[i][1], places=2)


    def test_run_zero_coverage(self):
        '''test run when all alleles have zero coverage'''
        mean_depth = 20
        error_rate = 0.01
        allele_combination_cov = {}
        allele_groups_dict = {'1': {0}, '2': {1}, '3': {0,1}, '4': {2}}
        allele_per_base_cov = [[0], [0,0]]
        gtyper = genotyper.Genotyper(mean_depth, error_rate, allele_combination_cov, allele_per_base_cov, allele_groups_dict)
        gtyper.run()
        self.assertEqual({'.'}, gtyper.genotype)
        self.assertEqual(0.0, gtyper.genotype_confidence)


