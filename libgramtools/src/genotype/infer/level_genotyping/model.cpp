#include "genotype/infer/level_genotyping/model.hpp"
#include "genotype/infer/allele_extracter.hpp"

using namespace gram::genotype::infer;
using namespace gram::genotype::infer::probabilities;

void LevelGenotyperModel::set_haploid_coverages(GroupedAlleleCounts const& gp_counts, AlleleId num_haplogroups){
    haploid_allele_coverages = PerAlleleCoverage(num_haplogroups, 0);
    singleton_allele_coverages = PerAlleleCoverage(num_haplogroups, 0);

    for (auto const& entry : gp_counts){
        for (auto const& allele_id : entry.first){
           haploid_allele_coverages.at(allele_id) += entry.second;
        }
        if (entry.first.size() == 1) {
            AlleleId id{entry.first.at(0)};
            singleton_allele_coverages.at(id) = entry.second;
        }
    }
}

void LevelGenotyperModel::assign_coverage_to_empty_alleles(allele_vector &alleles) {
    for (auto& allele : alleles){
        if (allele.sequence.size() == 0){
            auto assigned_cov = haploid_allele_coverages.at(allele.haplogroup);
            allele.pbCov = PerBaseCoverage{assigned_cov};
        }
    }
}

std::pair<double, double>
LevelGenotyperModel::compute_diploid_coverage(GroupedAlleleCounts const &gp_counts, AlleleIds ids,
                                              multiplicities const &haplogroup_multiplicities) {
    assert(ids.size() == 2);
    std::sort(ids.begin(), ids.end());

    if (computed_coverages.find(ids) != computed_coverages.end()){
        auto known_covs = computed_coverages.at(ids);
        return std::make_pair(known_covs.at(0), known_covs.at(1));
    }

    AlleleId first_allele_id = ids.at(0), second_allele_id = ids.at(1);

    double allele_1_cov = (double)(haploid_allele_coverages.at(first_allele_id));

    if (first_allele_id == second_allele_id) {
        // Case: two alleles from the same haplogroup are being compared
        if (haplogroup_multiplicities.at(first_allele_id)) {
            allele_1_cov /= 2;
            computed_coverages.insert({ids, allele_coverages{allele_1_cov, allele_1_cov}});
        }
        // Case: single homozygous coverage is evaluated
        else computed_coverages.insert({ids, allele_coverages{allele_1_cov}}); // Place homozygous cov once only

        return std::make_pair(allele_1_cov, allele_1_cov);
    }

    double allele_2_cov = (double)(haploid_allele_coverages.at(second_allele_id));
    bool has_first_allele, has_second_allele;
    CovCount shared_coverage{0};

    for (auto const& entry : gp_counts){
       has_first_allele = std::find(entry.first.begin(), entry.first.end(), first_allele_id) !=  entry.first.end();
       has_second_allele = std::find(entry.first.begin(), entry.first.end(), second_allele_id) !=  entry.first.end();

       if (has_first_allele && has_second_allele) shared_coverage += entry.second;
    }

    auto first_allele_specific_cov = allele_1_cov - shared_coverage,
        second_allele_specific_cov = allele_2_cov - shared_coverage;

    // Compute belonging factor (0 <= x <= 1)
    double allele1_belonging;
    if (first_allele_specific_cov == 0 && second_allele_specific_cov == 0) {
        // No allele specific coverage: equal allocation
        allele1_belonging = 0.5;
    }
    else allele1_belonging = first_allele_specific_cov / (first_allele_specific_cov + second_allele_specific_cov);

    allele_1_cov -= (1 - allele1_belonging) * shared_coverage;
    allele_2_cov -= allele1_belonging * shared_coverage;

    if (haplogroup_multiplicities.at(first_allele_id)) allele_1_cov /= 2;
    if (haplogroup_multiplicities.at(second_allele_id)) allele_2_cov /= 2;

    computed_coverages.insert({ids, allele_coverages{allele_1_cov, allele_2_cov}});
    return std::make_pair(allele_1_cov, allele_2_cov);
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

AlleleIds LevelGenotyperModel::get_haplogroups(allele_vector const &alleles, GtypedIndices const &gtype) const{
    AlleleIds result;
    for (auto const& index : gtype) {
        result.push_back(alleles.at(index).haplogroup);
    }
    std::sort(result.begin(), result.end());
    return result;
}


multiplicities LevelGenotyperModel::count_num_haplogroups(allele_vector const& alleles){
    std::map<AlleleId, std::size_t> haplo_counts;
    for (auto const& allele : alleles) {
        if (haplo_counts.find(allele.haplogroup) == haplo_counts.end())
            haplo_counts.insert({allele.haplogroup, 1});
        else haplo_counts.at(allele.haplogroup) += 1;
    }
    multiplicities m(haplo_counts.size(), false);
    for (auto const& entry : haplo_counts){
        if (entry.second > 1) m.at(entry.first) = true;
    }
    return m;
};

GtypedIndices LevelGenotyperModel::rescale_genotypes(GtypedIndices const& genotypes){
    // The rescaling does not affect allele 0, because that allele is always placed in the site
    std::unordered_map<GtypedIndex, GtypedIndex> rescaler{{0, 0}};

    GtypedIndices local_copy(genotypes), result;
    std::sort(local_copy.begin(), local_copy.end());

    GtypedIndex rescaling_index{1};

    for(auto const& gt : genotypes){
        if (rescaler.find(gt) == rescaler.end()) rescaler.insert(std::make_pair(gt, rescaling_index++));
        result.push_back(rescaler.at(gt));
    }

    return result;
}

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
        std::sort(combo_indices.begin(), combo_indices.end());
        result.push_back(combo_indices);
    } while (std::prev_permutation(v.begin(), v.end()));

    return result;
}


void LevelGenotyperModel::compute_haploid_log_likelihoods() {
    std::size_t allele_index{0};

    for ( auto const& allele : *alleles ){

       double cov_on_allele = haploid_allele_coverages.at(allele.haplogroup);
       auto cov_not_on_allele = total_coverage - cov_on_allele;

       double num_non_error_positions = count_credible_positions(l_stats->credible_cov_t, allele);
       double frac_non_error_positions = num_non_error_positions / allele.pbCov.size();

       double log_likelihood = (
              l_stats->poisson_full_depth.operator()(params{cov_on_allele}) +
              l_stats->log_mean_pb_error * cov_not_on_allele +
              frac_non_error_positions * l_stats->log_no_zero +
              (1 - frac_non_error_positions) * l_stats->mean_cov_depth * -1
              );

       likelihoods.insert({log_likelihood, GtypedIndices{allele_index}});
       ++allele_index;
    }
}

void LevelGenotyperModel::compute_homozygous_log_likelihoods(multiplicities const &haplogroup_multiplicities){
    std::size_t allele_index{0};

    for ( auto const& allele : *alleles ){

        auto coverages = compute_diploid_coverage(*gp_counts, AlleleIds{allele.haplogroup, allele.haplogroup},
                haplogroup_multiplicities);
        double cov_on_allele = coverages.first;
        auto cov_not_on_allele = total_coverage - cov_on_allele;

        cov_on_allele /= 2; // We will call poisson log-likelihood with half-depth on half the allele's coverage (twice).

        auto num_non_error_positions = count_credible_positions(l_stats->credible_cov_t, allele);
        double frac_non_error_positions = num_non_error_positions / allele.pbCov.size();

        double log_likelihood = (
                2 * l_stats->poisson_half_depth.operator()(params{cov_on_allele}) +
                l_stats->log_mean_pb_error * cov_not_on_allele +
                frac_non_error_positions * l_stats->log_no_zero +
                (1 - frac_non_error_positions) * l_stats->mean_cov_depth * -1
        );

        likelihoods.insert({log_likelihood, GtypedIndices{allele_index, allele_index}});
        ++allele_index;
    }
}

void LevelGenotyperModel::compute_heterozygous_log_likelihoods(multiplicities const &haplogroup_multiplicities) {
    GtypedIndices selected_indices;
    GtypedIndex index{0};
    for (auto const& allele : *alleles){
       if (singleton_allele_coverages.at(allele.haplogroup) != 0)
           selected_indices.push_back(index++);
    }

    if (selected_indices.size() < 2 ) return;

    auto all_diploid_combos = get_permutations(selected_indices, 2);

    for ( auto const& combo : all_diploid_combos ){
        Allele allele_1 = alleles->at(combo.at(0));
        Allele allele_2 = alleles->at(combo.at(1));
        AlleleIds haplogroups{
                allele_1.haplogroup,
                allele_2.haplogroup
        };
        auto coverages = compute_diploid_coverage(*gp_counts, haplogroups, haplogroup_multiplicities);
        auto allele_1_cov = coverages.first;
        auto allele_2_cov = coverages.second;

        auto allele_1_len = allele_1.pbCov.size();
        auto allele_2_len = allele_2.pbCov.size();

        double allele_1_frac_non_err_pos =
                count_credible_positions(l_stats->credible_cov_t, allele_1) /
                allele_1_len;

        double allele_2_frac_non_err_pos =
                count_credible_positions(l_stats->credible_cov_t, allele_2) /
                allele_2_len;

        double log_likelihood =
                l_stats->poisson_half_depth.operator()(params{allele_1_cov}) +
                l_stats->poisson_half_depth.operator()(params{allele_2_cov}) +
                (total_coverage - allele_1_cov - allele_2_cov) * l_stats->log_mean_pb_error +
                allele_1_frac_non_err_pos * l_stats->log_no_zero_half_depth +
                allele_2_frac_non_err_pos * l_stats->log_no_zero_half_depth +
                (1 - allele_1_frac_non_err_pos + 1 - allele_2_frac_non_err_pos) *
                        (-1 * l_stats->mean_cov_depth / 2);
        likelihoods.insert({log_likelihood, combo});
    }
}

LevelGenotyperModel::LevelGenotyperModel(allele_vector const *input_alleles, GroupedAlleleCounts const *gp_counts,
                                         Ploidy ploidy, likelihood_related_stats const *l_stats,
                                         bool ignore_ref_allele) :
gp_counts(gp_counts), ploidy(ploidy), l_stats(l_stats){
    assert(input_alleles->size() > 1);

    genotyped_site = std::make_shared<LevelGenotypedSite>();

    total_coverage = count_total_coverage(*gp_counts);
    genotyped_site->set_total_cov(total_coverage);
    if (total_coverage == 0 || l_stats->mean_cov_depth == 0){
        // Give the null genotyped site the first allele for extraction to be possible
        genotyped_site->set_alleles(allele_vector{input_alleles->at(0)});
        genotyped_site->make_null();
        return;
    }

    allele_vector used_alleles;
    if (ignore_ref_allele)
        used_alleles = allele_vector{input_alleles->begin() + 1, input_alleles->end()};
    else used_alleles = *input_alleles;
    this->alleles = &used_alleles;

    auto haplogroup_multiplicities = count_num_haplogroups(*alleles);
    genotyped_site->set_num_haplogroups(haplogroup_multiplicities.size());
    set_haploid_coverages(*gp_counts, haplogroup_multiplicities.size());
    assign_coverage_to_empty_alleles(*alleles);

    if (ploidy == Ploidy::Haploid) compute_haploid_log_likelihoods();
    else if (ploidy == Ploidy::Diploid){
        compute_homozygous_log_likelihoods(haplogroup_multiplicities);
        compute_heterozygous_log_likelihoods(haplogroup_multiplicities);
    }

    auto it = likelihoods.begin();
    auto chosen_gt = it->second;
    auto best_likelihood = it->first;
    ++it;
    auto gt_confidence = best_likelihood - it->first;

    // Correct the genotype indices if needed
    if (ignore_ref_allele) for (auto& gt : chosen_gt) gt++;

    auto chosen_alleles = genotyped_site->get_unique_genotyped_alleles(*input_alleles, chosen_gt);
    // Get haplotypes and coverages
    auto chosen_haplotypes = get_haplogroups(*input_alleles, chosen_gt);
    allele_coverages allele_covs;
    if (ploidy == Ploidy::Haploid) allele_covs =
            allele_coverages{(double)haploid_allele_coverages.at(chosen_haplotypes.at(0))};
    else if (ploidy == Ploidy::Diploid) allele_covs = computed_coverages.at(chosen_haplotypes);

    auto rescaled_gt = rescale_genotypes(chosen_gt);
    // Add the REF allele to the set of chosen alleles if it was not called
    if (rescaled_gt.at(0) != 0) chosen_alleles = prepend(chosen_alleles, input_alleles->at(0));
    if (chosen_haplotypes.at(0) != 0) {
        auto ref_cov = (double)singleton_allele_coverages.at(0);
        allele_covs = prepend(allele_covs, ref_cov);
    }

    genotyped_site->set_alleles(chosen_alleles);
    genotyped_site->set_genotype(rescaled_gt, gt_confidence);
}

