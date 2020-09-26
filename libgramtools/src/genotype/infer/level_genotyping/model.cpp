#include "genotype/infer/level_genotyping/model.hpp"

#include "genotype/infer/allele_extracter.hpp"

using namespace gram::genotype::infer;
using namespace gram::genotype::infer::probabilities;

std::optional<Allele> check_for_duplicates(allele_vector const& input_alleles) {
  std::set<Allele> nodups;
  for (auto const& allele : input_alleles) {
    if (nodups.find(allele) == nodups.end())
      nodups.insert(allele);
    else
      return allele;
  }
  return std::nullopt;
}

LevelGenotyperModel::LevelGenotyperModel(ModelData& input_data)
    : data(input_data) {
  assert(data.input_alleles.size() > 1);
  auto const& ref_allele = data.input_alleles.at(0);
  genotyped_site = std::make_shared<LevelGenotypedSite>();

  auto haplogroup_multiplicities =
      get_haplogroup_multiplicities(data.input_alleles);
  // Used in invalidation process, required set early
  genotyped_site->set_num_haplogroups(haplogroup_multiplicities.size());

  auto const duplicate_allele = check_for_duplicates(data.input_alleles);
  if (duplicate_allele) genotyped_site->set_filter("AMBIG");

  total_coverage = count_total_coverage(data.gp_counts);
  if (total_coverage == 0 || data.l_stats->data_params.mean_cov == 0) {
    // Null gt site still gets ref allele, for reporting and for allele
    // extraction to work
    genotyped_site->set_alleles(allele_vector{ref_allele});
    genotyped_site->make_null();
    return;
  }

  set_haploid_coverages(data.gp_counts, haplogroup_multiplicities.size());

  allele_vector used_alleles(data.input_alleles);
  assign_coverage_to_empty_alleles(used_alleles);

  if (data.ploidy == Ploidy::Haploid)
    compute_haploid_log_likelihoods(used_alleles);
  else if (data.ploidy == Ploidy::Diploid) {
    compute_homozygous_log_likelihoods(used_alleles, haplogroup_multiplicities);
    compute_heterozygous_log_likelihoods(used_alleles,
                                         haplogroup_multiplicities);
  }

  CallGenotype(data.input_alleles, haplogroup_multiplicities, data.ploidy);
}

void LevelGenotyperModel::set_haploid_coverages(
    GroupedAlleleCounts const& input_gp_counts, AlleleId num_haplogroups) {
  haploid_allele_coverages = PerAlleleCoverage(num_haplogroups, 0);
  singleton_allele_coverages = PerAlleleCoverage(num_haplogroups, 0);

  for (auto const& entry : input_gp_counts) {
    for (auto const& allele_id : entry.first) {
      haploid_allele_coverages.at(allele_id) += entry.second;
    }
    if (entry.first.size() == 1) {
      AlleleId id{entry.first.at(0)};
      singleton_allele_coverages.at(id) = entry.second;
    }
  }
}

void LevelGenotyperModel::assign_coverage_to_empty_alleles(
    allele_vector& input_alleles) {
  for (auto& allele : input_alleles) {
    if (allele.sequence.empty()) {
      auto assigned_cov = haploid_allele_coverages.at(allele.haplogroup);
      allele.pbCov = PerBaseCoverage{assigned_cov};
    }
  }
}

CovPair LevelGenotyperModel::diploid_cov_same_haplogroup(
    AlleleIds const& haplogroups) {
  auto hapg = haplogroups.at(0);
  auto allele_cov = (double)(haploid_allele_coverages.at(hapg));
  allele_cov /= 2;
  computed_coverages.insert(
      {haplogroups, allele_coverages{allele_cov, allele_cov}});
  return std::make_pair(allele_cov, allele_cov);
}

std::pair<double, double> LevelGenotyperModel::diploid_cov_different_haplogroup(
    GroupedAlleleCounts const& gp_counts, AlleleIds const& ids,
    multiplicities const& hap_mults) {
  AlleleId allele_1_id = ids.at(0), allele_2_id = ids.at(1);
  auto allele_1_cov = (double)(haploid_allele_coverages.at(allele_1_id));
  auto allele_2_cov = (double)(haploid_allele_coverages.at(allele_2_id));
  bool has_first_allele, has_second_allele;
  CovCount shared_coverage{0};

  for (auto const& entry : gp_counts) {
    has_first_allele = std::find(entry.first.begin(), entry.first.end(),
                                 allele_1_id) != entry.first.end();
    has_second_allele = std::find(entry.first.begin(), entry.first.end(),
                                  allele_2_id) != entry.first.end();
    if (has_first_allele && has_second_allele) shared_coverage += entry.second;
  }

  auto first_allele_specific_cov = allele_1_cov - shared_coverage,
       second_allele_specific_cov = allele_2_cov - shared_coverage;

  // Compute belonging factor (0 <= x <= 1)
  double allele1_belonging;
  if (first_allele_specific_cov == 0 && second_allele_specific_cov == 0) {
    // No allele specific coverage: equal allocation
    allele1_belonging = 0.5;
  } else
    allele1_belonging =
        first_allele_specific_cov /
        (first_allele_specific_cov + second_allele_specific_cov);

  allele_1_cov -= (1 - allele1_belonging) * shared_coverage;
  allele_2_cov -= allele1_belonging * shared_coverage;

  if (hap_mults.at(allele_1_id)) allele_1_cov /= 2;
  if (hap_mults.at(allele_2_id)) allele_2_cov /= 2;

  computed_coverages.insert(
      {ids, allele_coverages{allele_1_cov, allele_2_cov}});
  return std::make_pair(allele_1_cov, allele_2_cov);
}

std::pair<double, double> LevelGenotyperModel::compute_diploid_coverage(
    GroupedAlleleCounts const& gp_counts, AlleleIds haplogroups,
    multiplicities const& haplogroup_multiplicities) {
  assert(haplogroups.size() == 2);
  // Below line so that `computed_coverages` map lookup/insertion is
  // insensitive to order.
  std::sort(haplogroups.begin(), haplogroups.end());

  if (computed_coverages.find(haplogroups) != computed_coverages.end()) {
    auto known_covs = computed_coverages.at(haplogroups);
    return std::make_pair(known_covs.at(0), known_covs.at(1));
  }

  if (haplogroups.at(0) == haplogroups.at(1))
    return diploid_cov_same_haplogroup(haplogroups);
  else
    return diploid_cov_different_haplogroup(gp_counts, haplogroups,
                                            haplogroup_multiplicities);
}

double LevelGenotyperModel::fraction_noncredible_positions(
    Allele const& allele) {
  auto& credible_cov_threshold = data.l_stats->credible_cov_t;
  double num_noncredible_pos{0.};
  for (auto const& pb_cov : allele.pbCov) {
    if (pb_cov < credible_cov_threshold) ++num_noncredible_pos;
  }
  double result = num_noncredible_pos / allele.pbCov.size();
  return result;
}

std::size_t LevelGenotyperModel::count_total_coverage(
    GroupedAlleleCounts const& gp_counts) {
  std::size_t total_cov{0};
  for (auto const& entry : gp_counts) total_cov += entry.second;
  return total_cov;
}

AlleleIds LevelGenotyperModel::get_haplogroups(
    allele_vector const& alleles, GtypedIndices const& gtype) const {
  AlleleIds result;
  for (auto const& index : gtype) {
    result.push_back(alleles.at(index).haplogroup);
  }
  std::sort(result.begin(), result.end());
  return result;
}

multiplicities LevelGenotyperModel::get_haplogroup_multiplicities(
    allele_vector const& input_alleles) {
  std::map<AlleleId, std::size_t> haplo_counts;
  for (auto const& allele : input_alleles) {
    if (haplo_counts.find(allele.haplogroup) == haplo_counts.end())
      haplo_counts.insert({allele.haplogroup, 1});
    else
      haplo_counts.at(allele.haplogroup) += 1;
  }
  multiplicities m(haplo_counts.size(), false);
  for (auto const& entry : haplo_counts) {
    if (entry.second > 1) m.at(entry.first) = true;
  }
  return m;
}

GtypedIndices LevelGenotyperModel::rescale_genotypes(
    GtypedIndices const& genotypes) {
  // The rescaling does not affect allele 0, because that allele is always
  // placed in the site
  std::unordered_map<GtypedIndex, GtypedIndex> rescaler{{0, 0}};

  GtypedIndices local_copy(genotypes), result;
  std::sort(local_copy.begin(), local_copy.end());

  GtypedIndex rescaling_index{1};

  for (auto const& gt : genotypes) {
    if (rescaler.find(gt) == rescaler.end())
      rescaler.insert(std::make_pair(gt, rescaling_index++));
    result.push_back(rescaler.at(gt));
  }
  return result;
}

std::vector<GtypedIndices> LevelGenotyperModel::get_permutations(
    const GtypedIndices& indices, std::size_t const subset_size) {
  // Credit: https://stackoverflow.com/a/9430993/12519542
  auto num_elements = indices.size();
  if (subset_size > num_elements) return {};

  std::vector<bool> v(num_elements);
  std::fill(v.begin(), v.begin() + subset_size, true);

  std::vector<GtypedIndices> result;
  do {
    GtypedIndices combo_indices;
    for (int i = 0; i < num_elements; ++i) {
      if (v[i]) combo_indices.push_back(indices.at(i));
    }
    std::sort(combo_indices.begin(), combo_indices.end());
    result.push_back(combo_indices);
  } while (std::prev_permutation(v.begin(), v.end()));

  return result;
}

void LevelGenotyperModel::add_likelihood(allele_vector const& alleles,
                                         double const& incompatible_coverage,
                                         GtypedIndices const& allele_indices) {
  double log_likelihood =
      incompatible_coverage * data.l_stats->log_mean_pb_error;
  uint8_t stop_index;
  switch (data.ploidy) {
    case Ploidy::Haploid:
      stop_index = 1;
      break;
    case Ploidy::Diploid:
      stop_index = 2;
      break;
    default:
      throw UnsupportedPloidy("");
  }

  assert(alleles.size() == stop_index);
  assert(allele_indices.size() == stop_index);
  double compatible_coverage, gap_penalty;
  for (uint8_t i = 0; i < stop_index; i++) {
    auto const& allele = alleles.at(i);
    compatible_coverage = allele.get_average_cov();
    gap_penalty = fraction_noncredible_positions(allele);
    log_likelihood +=
        data.l_stats->pmf_full_depth->operator()(params{compatible_coverage});
    log_likelihood += gap_penalty * data.l_stats->log_zero;
  }

  likelihoods.insert({log_likelihood, allele_indices});
}

void LevelGenotyperModel::compute_haploid_log_likelihoods(
    allele_vector const& input_alleles) {
  GtypedIndex allele_index{-1};

  for (auto const& allele : input_alleles) {
    ++allele_index;
    if (allele_index == 0 && ignore_ref_allele()) continue;
    auto haploid_cov = haploid_allele_coverages.at(allele.haplogroup);
    auto incompatible_coverage = total_coverage - haploid_cov;
    add_likelihood(allele_vector{allele}, incompatible_coverage,
                   GtypedIndices{allele_index});
  }
}

void LevelGenotyperModel::compute_homozygous_log_likelihoods(
    allele_vector const& input_alleles,
    multiplicities const& haplogroup_multiplicities) {
  GtypedIndex allele_index{-1};

  for (auto const& allele : input_alleles) {
    ++allele_index;
    if (allele_index == 0 && ignore_ref_allele()) continue;
    auto coverages = compute_diploid_coverage(
        data.gp_counts, AlleleIds{allele.haplogroup, allele.haplogroup},
        haplogroup_multiplicities);
    auto incompatible_coverage =
        total_coverage - coverages.first - coverages.second;

    add_likelihood(allele_vector{allele, allele}, incompatible_coverage,
                   GtypedIndices{allele_index, allele_index});
  }
}

void LevelGenotyperModel::compute_heterozygous_log_likelihoods(
    allele_vector const& input_alleles,
    multiplicities const& haplogroup_multiplicities) {
  GtypedIndices selected_indices;
  GtypedIndex allele_index{-1};
  for (auto const& allele : input_alleles) {
    ++allele_index;
    if (allele_index == 0 && ignore_ref_allele()) continue;
    if (singleton_allele_coverages.at(allele.haplogroup) != 0)
      selected_indices.push_back(allele_index);
  }

  if (selected_indices.size() < 2) return;

  auto all_diploid_combos = get_permutations(selected_indices, 2);

  for (auto const& combo : all_diploid_combos) {
    Allele allele_1 = input_alleles.at(combo.at(0));
    Allele allele_2 = input_alleles.at(combo.at(1));
    AlleleIds haplogroups{allele_1.haplogroup, allele_2.haplogroup};
    auto coverages = compute_diploid_coverage(data.gp_counts, haplogroups,
                                              haplogroup_multiplicities);
    auto incompatible_coverage =
        total_coverage - coverages.first - coverages.second;

    add_likelihood(allele_vector{allele_1, allele_2}, incompatible_coverage,
                   combo);
  }
}

void LevelGenotyperModel::add_next_best_alleles(
    allele_vector const& input_alleles, GtypedIndices const& chosen_gt,
    GtypedIndices const& next_best_gt) {
  auto& chosen_allele_for_cov = input_alleles.at(chosen_gt.at(0));
  auto& next_best_allele_for_cov = input_alleles.at(next_best_gt.at(0));

  bool low_total_cov = total_coverage < data.l_stats->data_params.mean_cov / 4;
  bool low_relative_cov =
      haploid_allele_coverages.at(chosen_allele_for_cov.haplogroup) <
      haploid_allele_coverages.at(next_best_allele_for_cov.haplogroup) * 2;

  if (low_total_cov || low_relative_cov) {
    std::set<GtypedIndex> next_best{next_best_gt.begin(), next_best_gt.end()};
    for (auto const& gt : chosen_gt) {
      // Can need erasures in diploid genotyping
      if (next_best.find(gt) != next_best.end()) next_best.erase(gt);
    }
    allele_vector result;
    for (auto const& gt : next_best) {
      auto added_allele = input_alleles.at(gt);
      added_allele.nesting_consistent = false;
      result.push_back(added_allele);
    }
    genotyped_site->set_extra_alleles(result);
  }
}

void LevelGenotyperModel::add_all_best_alleles(
    allele_vector const& input_alleles, GtypedIndices const& chosen_gt,
    GtypedIndices const& next_best_gt) {
  std::set<GtypedIndex> all_best{next_best_gt.begin(), next_best_gt.end()};
  all_best.insert(chosen_gt.begin(), chosen_gt.end());
  allele_vector result;
  for (auto const& gt : all_best) result.push_back(input_alleles.at(gt));
  genotyped_site->set_extra_alleles(result);
}

likelihood_map::const_iterator LevelGenotyperModel::ChooseMaxLikelihood(
    likelihood_map const& likelihoods, allele_vector const& alleles) {
  if (likelihoods.size() < 2)
    throw IncorrectGenotyping(
        "Less than 2 alleles have a likelihood.\n"
        "Allele extraction bug?");
  auto it = likelihoods.begin();
  while (it != likelihoods.end()) {
    bool inconsistent = false;
    auto chosen_gt = it->second;
    for (auto const& gt : chosen_gt) {
      if (!alleles.at(gt).nesting_consistent) {
        inconsistent = true;
        break;
      }
    }
    if (inconsistent)
      ++it;
    else
      break;
  }
  if (std::distance(it, likelihoods.end()) < 2)
    throw IncorrectGenotyping(
        "Fewer than 2 alleles are consistent with child"
        "sites.\nAllele extraction bug?");

  return it;
}

void LevelGenotyperModel::CallGenotype(allele_vector const& input_alleles,
                                       multiplicities hap_mults,
                                       Ploidy const ploidy) {
  auto const& ref_allele = input_alleles.at(0);
  auto it = ChooseMaxLikelihood(likelihoods, input_alleles);
  auto best_likelihood = it->first;
  auto chosen_gt = it->second;
  ++it;
  auto gt_confidence = best_likelihood - it->first;
  auto next_best_gt = it->second;

  if (gt_confidence == 0.) {
    genotyped_site->set_alleles(allele_vector{ref_allele});
    genotyped_site->make_null();
    add_all_best_alleles(input_alleles, chosen_gt, next_best_gt);
    return;
  } else
    add_next_best_alleles(input_alleles, chosen_gt, next_best_gt);

  auto chosen_alleles =
      genotyped_site->get_unique_genotyped_alleles(input_alleles, chosen_gt);
  // Get haplotypes and coverages
  auto chosen_haplotypes = get_haplogroups(input_alleles, chosen_gt);
  allele_coverages allele_covs;
  if (ploidy == Ploidy::Haploid)
    allele_covs = allele_coverages{
        (double)haploid_allele_coverages.at(chosen_haplotypes.at(0))};
  else if (ploidy == Ploidy::Diploid) {
    allele_covs = computed_coverages.at(chosen_haplotypes);
    // If homozygous call, give the single allele all the coverage
    if (chosen_gt.at(0) == chosen_gt.at(1))
      allele_covs = allele_coverages{allele_covs.at(0) + allele_covs.at(1)};
  }

  // Because we only output the called alleles (+ REF), re-express genotype
  // indices to those.
  auto rescaled_gt = rescale_genotypes(chosen_gt);

  // Add the REF allele (and its coverage) to the set of chosen alleles if it
  // was not called
  if (rescaled_gt.at(0) != 0) {
    chosen_alleles = prepend(chosen_alleles, ref_allele);
    auto ref_cov = (double)singleton_allele_coverages.at(0);
    if (hap_mults.at(0)) ref_cov /= 2;
    allele_covs = prepend(allele_covs, ref_cov);
  }

  assert(chosen_alleles.size() == allele_covs.size());
  genotyped_site->populate_site(gtype_information{
      chosen_alleles, rescaled_gt, allele_covs, total_coverage,
      genotyped_site->get_genotyped_haplogroups(chosen_alleles, rescaled_gt)});
  genotyped_site->set_gt_conf(gt_confidence);

  if (data.debug) {
    std::string debug_info{"\tnext_best_seq: "};
    for (auto const& gt : next_best_gt) {
      debug_info.append(input_alleles.at(gt).sequence);
      debug_info.append(",");
    }
    debug_info.append("\tnext_best_cov: ");
    auto next_best_haplotypes = get_haplogroups(input_alleles, next_best_gt);
    for (auto const& hapg : next_best_haplotypes) {
      debug_info.append(std::to_string(haploid_allele_coverages.at(hapg)));
      debug_info.append(",");
    }
    genotyped_site->set_debug_info(debug_info);
  }
}

LevelGenotyperModel::LevelGenotyperModel(
    likelihood_related_stats const& input_l_stats,
    PerAlleleCoverage const& input_covs,
    likelihood_map const& input_likelihoods)
    : likelihoods(input_likelihoods) {
  data.l_stats = &input_l_stats;
  genotyped_site = std::make_shared<LevelGenotypedSite>();

  haploid_allele_coverages = input_covs;
  singleton_allele_coverages = input_covs;
  total_coverage = 0;
  for (auto const& entry : input_covs) total_coverage += entry;
}
