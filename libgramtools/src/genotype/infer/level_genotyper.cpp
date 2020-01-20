#include "genotype/infer/genotyping_models.hpp"

using namespace gram::genotype::infer;
using namespace gram::genotype::infer::probabilities;

void LevelGenotyperModel::set_haploid_coverages(GroupedAlleleCounts const& gp_counts, AlleleId num_haplogroups){
    haploid_allele_coverages = PerAlleleCoverage(num_haplogroups, 0);

    for (auto const& entry : gp_counts){
        for (auto const& allele_id : entry.first){
           haploid_allele_coverages.at(allele_id) += entry.second;
        }
        if (entry.first.size() == 1) {
            AlleleId id{entry.first.at(0)};
            singleton_allele_coverages.insert(id);
        }
    }
}

std::pair<float, float> LevelGenotyperModel::compute_diploid_coverage(GroupedAlleleCounts const& gp_counts, AlleleIds ids){
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


numCredibleCounts LevelGenotyperModel::count_credible_positions(CovCount const& credible_cov_t, Allele const& allele){
    numCredibleCounts c{0};
    for (auto const& pb_cov : allele.pbCov){
        if (pb_cov >= credible_cov_t) ++c;
    }
    return c;
}

std::size_t LevelGenotyperModel::count_total_coverage(GroupedAlleleCounts const& gp_counts){
    std::size_t total_cov{0};
    for (auto const& entry : gp_counts) total_cov += entry.second;
    return total_cov;
}

std::size_t LevelGenotyperModel::count_num_haplogroups(allele_vector const& alleles){
    std::set<AlleleId> unique_haplos;
    for (auto const& allele : alleles) {
        if (unique_haplos.find(allele.haplogroup) == unique_haplos.end())
            unique_haplos.insert(allele.haplogroup);
    }
    return unique_haplos.size();
};

std::vector<GtypedIndices> LevelGenotyperModel::get_permutations(const GtypedIndices &indices, std::size_t const subset_size){
    auto num_elements = indices.size();
    if (subset_size > num_elements) return {};

    std::vector<bool> v(num_elements);
    std::fill(v.begin(), v.begin() + subset_size, true);

    std::vector<GtypedIndices> result;
    do {
        GtypedIndices combo_indices;
        for (int i = 0; i < num_elements; ++i){
            if (v[i]) combo_indices.push_back(indices.at(i));
        }
        result.push_back(combo_indices);
    } while (std::prev_permutation(v.begin(), v.end()));

    return result;
}


void LevelGenotyperModel::compute_homogeneous_log_likelihoods(bool haploid) {
    std::size_t allele_index{0};

    for ( auto const& allele : *alleles ){

       double cov_on_allele = haploid_allele_coverages.at(allele.haplogroup);
       auto cov_not_on_allele = total_coverage - cov_on_allele;

       auto num_non_error_positions = count_credible_positions(l_stats->credible_cov_t, allele);
       double frac_non_error_positions = num_non_error_positions / allele.pbCov.size();

       double likelihood = (
              l_stats->poisson_full_depth.operator()(params{cov_on_allele}) +
              l_stats->log_mean_pb_error * cov_not_on_allele +
              frac_non_error_positions * l_stats->log_no_zero +
              (1 - frac_non_error_positions) * l_stats->mean_cov_depth * -1
              );

       if (haploid) likelihoods.insert({likelihood, GtypedIndices{allele_index}});
       else likelihoods.insert({likelihood, GtypedIndices{allele_index, allele_index}});

       ++allele_index;
    }
}

void LevelGenotyperModel::compute_heterozygous_log_likelihoods(){
    GtypedIndices selected_indices;
    GtypedIndex index{0};
    for (auto const& allele : *alleles){
       if (singleton_allele_coverages.find(allele.haplogroup) != singleton_allele_coverages.end())
           selected_indices.push_back(index++);
    }

    if (selected_indices.size() < 2 ) return;

    auto all_diploid_combos = get_permutations(selected_indices, 2);

    for ( auto const& combo : all_diploid_combos ){
        Allele first_allele = alleles->at(combo.at(0));
        Allele second_allele = alleles->at(combo.at(1));
        AlleleIds haplogroups{
                first_allele.haplogroup,
                second_allele.haplogroup
        };
       auto coverages = compute_diploid_coverage(*gp_counts, haplogroups);
    }
}

LevelGenotyperModel::LevelGenotyperModel(allele_vector const *alleles, GroupedAlleleCounts const *gp_counts,
                                         Ploidy ploidy,
                                         likelihood_related_stats const *l_stats) :
alleles(alleles), gp_counts(gp_counts), ploidy(ploidy), l_stats(l_stats){
    genotyped_site = std::make_shared<LevelGenotypedSite>();

    total_coverage = count_total_coverage(*gp_counts);
    if (total_coverage == 0 || l_stats->mean_cov_depth == 0){
        genotyped_site->make_null();
        return;
    }

    auto num_haplogroups = count_num_haplogroups(*alleles);
    set_haploid_coverages(*gp_counts, num_haplogroups);

    if (ploidy == Ploidy::Haploid) compute_homogeneous_log_likelihoods(true);
    else if (ploidy == Ploidy::Diploid){
       compute_homogeneous_log_likelihoods(false);
       compute_heterozygous_log_likelihoods();
    }

    auto it = likelihoods.begin();
    auto chosen_gt = it->second;
    auto best_likelihood = it->first;
    ++it;
    auto gt_confidence = best_likelihood - it->first;

    genotyped_site->set_genotype(chosen_gt, gt_confidence);
}

