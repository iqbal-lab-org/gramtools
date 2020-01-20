#include "genotype/infer/genotyping_models.hpp"

using namespace gram::genotype::infer;
using namespace gram::genotype::infer::probabilities;

void LevelGenotyper::set_haploid_coverages(GroupedAlleleCounts const& gp_counts, AlleleId num_haplogroups){
    haploid_allele_coverages = singleton_allele_coverages = PerAlleleCoverage(num_haplogroups, 0);

    for (auto const& entry : gp_counts){
        for (auto const& allele_id : entry.first){
           haploid_allele_coverages.at(allele_id) += entry.second;
        }
        if (entry.first.size() == 1) {
            AlleleId id{entry.first.at(0)};
            singleton_allele_coverages.at(id) += entry.second;
        }
    }
}

std::pair<float, float> LevelGenotyper::compute_diploid_coverage(GroupedAlleleCounts const& gp_counts, AlleleIds ids){
    assert(ids.size() == 2);
    AlleleId first_allele_id = ids.at(0), second_allele_id = ids.at(1);
    float first_allele_coverage = (float)(haploid_allele_coverages.at(first_allele_id));
    float second_allele_coverage = (float)(haploid_allele_coverages.at(second_allele_id));

    bool has_first_allele, has_second_allele;
    CovCount shared_coverage{0};
    for (auto const& entry : gp_counts){
       has_first_allele = std::find(entry.first.begin(), entry.first.end(), first_allele_id) !=  entry.first.end();
       has_second_allele = std::find(entry.first.begin(), entry.first.end(), second_allele_id) !=  entry.first.end();

       if (has_first_allele && has_second_allele) shared_coverage += entry.second;
    }

    auto first_allele_specific_cov = first_allele_coverage - shared_coverage,
        second_allele_specific_cov = second_allele_coverage - shared_coverage;

    // Compute belonging factor (0 <= x <= 1)
    float allele1_belonging;
    if (first_allele_specific_cov == 0 && second_allele_specific_cov == 0) {
        // No allele specific coverage: equal allocation
        allele1_belonging = 0.5;
    }
    else allele1_belonging = first_allele_specific_cov / (first_allele_specific_cov + second_allele_specific_cov);

    first_allele_coverage -= (1 - allele1_belonging) * shared_coverage;
    second_allele_coverage -= allele1_belonging * shared_coverage;

    return std::make_pair(first_allele_coverage, second_allele_coverage);
}


numCredibleCounts LevelGenotyper::count_credible_positions(CovCount const& credible_cov_t, Allele const& allele){
    numCredibleCounts c{0};
    for (auto const& pb_cov : allele.pbCov){
        if (pb_cov >= credible_cov_t) ++c;
    }
    return c;
}

LevelGenotyper::LevelGenotyper(allele_vector const* alleles, GroupedAlleleCounts const* gp_counts, Ploidy ploidy,
poisson_pmf_ptr poisson_prob, likelihood_related_stats const* l_stats) :
alleles(alleles), gp_counts(gp_counts), ploidy(ploidy), poisson_prob(poisson_prob), l_stats(l_stats){
    ;
}
