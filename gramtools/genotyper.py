import itertools
import math
import operator

from scipy.stats import poisson


class Error (Exception): pass

class Genotyper:
    def __init__(self, mean_depth, error_rate, allele_combination_cov, allele_per_base_cov, allele_groups_dict):
        self.mean_depth = mean_depth
        self.error_rate = error_rate
        self.allele_combination_cov = allele_combination_cov
        self.allele_per_base_cov = allele_per_base_cov
        self.allele_groups_dict = allele_groups_dict
        self.likelihoods = None
        self.genotype = None
        self.genotype_confidence = None
        self.singleton_alleles_cov = {}


    @classmethod
    def _singleton_alleles_and_coverage(cls, allele_combination_cov, allele_groups_dict):
        singleton_alleles = {}
        for allele_key in allele_combination_cov:
            if len(allele_groups_dict[allele_key]) == 1:
                allele = allele_groups_dict[allele_key].pop()
                allele_groups_dict[allele_key].add(allele)
                singleton_alleles[allele] = allele_combination_cov[allele_key]
        return singleton_alleles


    @classmethod
    def _total_coverage(cls, allele_combination_cov):
        return sum(allele_combination_cov.values())


    @classmethod
    def _coverage_of_one_haploid_allele(cls, allele, allele_combination_cov, allele_groups_dict):
        coverage = 0
        for allele_key in allele_combination_cov:
            if allele in allele_groups_dict[allele_key]:
                coverage += allele_combination_cov[allele_key]
        return coverage


    @classmethod
    def _coverage_of_diploid_alleles(cls, allele1, allele2, allele_combination_cov, allele_groups_dict, singleton_alleles_cov):
        assert allele1 in singleton_alleles_cov and allele2 in singleton_alleles_cov
        allele1_cov = singleton_alleles_cov[allele1]
        allele2_cov = singleton_alleles_cov[allele2]
        allele1_total_cov = 0
        allele2_total_cov = 0

        for allele_key in allele_combination_cov:
            allele_combination = allele_groups_dict[allele_key]
            if allele_combination.issuperset({allele1, allele2}):
                allele1_total_cov += (allele1_cov / (allele1_cov + allele2_cov) ) * allele_combination_cov[allele_key]
                allele2_total_cov += (allele2_cov / (allele1_cov + allele2_cov) ) * allele_combination_cov[allele_key]
            elif allele1 in allele_combination:
                assert allele2 not in allele_combination
                allele1_total_cov += allele_combination_cov[allele_key]
            elif allele2 in allele_combination:
                assert allele1 not in allele_combination
                allele2_total_cov += allele_combination_cov[allele_key]

        return allele1_total_cov, allele2_total_cov


    @classmethod
    def _log_likelihood_homozygous(cls, mean_depth, allele_depth, total_depth, error_rate, allele_length, non_zeros):
        return sum([
            -mean_depth * (1 + allele_length - non_zeros),
            allele_depth * math.log(mean_depth),
            -math.lgamma(allele_depth + 1),
            (total_depth - allele_depth) * math.log(error_rate),
            non_zeros * math.log(1 - poisson.pmf(0, mean_depth)),
        ])



    @classmethod
    def _log_likelihood_heterozygous(cls, mean_depth, allele_depth1, allele_depth2, total_depth,
            error_rate, allele_length1, allele_length2, non_zeros1, non_zeros2):
        return sum([
            -mean_depth * (1 + 0.5 * (allele_length1 + allele_length2 - non_zeros1 - non_zeros2)),
            (allele_depth1 + allele_depth2) * math.log(0.5 * mean_depth),
            -math.lgamma(allele_depth1 + 1),
            -math.lgamma(allele_depth2 + 1),
            (total_depth - allele_depth1 - allele_depth2) * math.log(error_rate),
            (non_zeros1 + non_zeros2) * math.log(1 - poisson.pmf(0, 0.5 * mean_depth)),
        ])


    def _calculate_log_likelihoods(self):
        '''Makes a list of tuples: ( (allele(s) tuple), log likelihood).
        List is sorted from most to least likely'''
        self.likelihoods = []
        total_depth = sum(self.allele_combination_cov.values())

        for allele_number, per_base_cov in enumerate(self.allele_per_base_cov):
            allele_depth = Genotyper._coverage_of_one_haploid_allele(allele_number, self.allele_combination_cov, self.allele_groups_dict)
            allele_length = len(per_base_cov)
            non_zeros = len(per_base_cov) - per_base_cov.count(0)

            log_likelihood = Genotyper._log_likelihood_homozygous(
                    self.mean_depth,
                    allele_depth,
                    total_depth,
                    self.error_rate,
                    allele_length,
                    non_zeros,
            )
            self.likelihoods.append(({allele_number, allele_number}, log_likelihood))

        self.singleton_alleles_cov = Genotyper._singleton_alleles_and_coverage(self.allele_combination_cov, self.allele_groups_dict)

        for (allele_number1, allele_number2) in itertools.combinations(self.singleton_alleles_cov.keys(), 2):
            allele1_depth = self.singleton_alleles_cov[allele_number1]
            allele2_depth = self.singleton_alleles_cov[allele_number2]
            allele1_length = len(self.allele_per_base_cov[allele_number1])
            allele2_length = len(self.allele_per_base_cov[allele_number2])
            non_zeros1 = allele1_length - self.allele_per_base_cov[allele_number1].count(0)
            non_zeros2 = allele2_length - self.allele_per_base_cov[allele_number2].count(0)

            log_likelihood = Genotyper._log_likelihood_heterozygous(
                    self.mean_depth,
                    allele1_depth,
                    allele2_depth,
                    total_depth,
                    self.error_rate,
                    allele1_length,
                    allele2_length,
                    non_zeros1,
                    non_zeros2,
            )
            self.likelihoods.append((set([allele_number1, allele_number2]), log_likelihood))

        self.likelihoods.sort(key=operator.itemgetter(1),reverse=True)


    def run(self):
        if len(self.allele_combination_cov) == 0 or Genotyper._total_coverage(self.allele_combination_cov) == 0:
            self.genotype = {'.'}
            self.genotype_confidence = 0.0
        else:
            self._calculate_log_likelihoods()
            assert self.likelihoods is not None and len(self.likelihoods) > 1
            self.genotype, best_log_likelihood = self.likelihoods[0]
            self.genotype_confidence = round(best_log_likelihood - self.likelihoods[1][1], 2)

