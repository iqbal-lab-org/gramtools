#include "genotype/infer/level_genotyping/site.hpp"
#include "genotype/infer/output_specs/coverage_common.hpp"

using namespace gram::genotype::infer;

gram::AlleleIds const LevelGenotypedSite::get_nonGenotyped_haplogroups() const {
  assert(gtype_info.alleles.size() > 0);
  assert(num_haplogroups > 0);
  AlleleIds result;

  AlleleIdSet genotyped_haplogroups;
  if (!is_null())
    for (auto const& gt : gtype_info.genotype)
      genotyped_haplogroups.insert(gtype_info.alleles.at(gt).haplogroup);

  for (AlleleId i{0}; i < num_haplogroups; i++) {
    if (genotyped_haplogroups.find(i) == genotyped_haplogroups.end())
      result.push_back(i);
  }
  return result;
}

site_entries LevelGenotypedSite::get_model_specific_entries() {
  site_entry<double> gt_conf_entry{"FORMAT", "GT_CONF", {gt_conf}, true};
  site_entry<double> gt_conf_percentile_entry{
      "FORMAT", "GT_CONF_PERCENTILE", {gt_conf_percentile}, true};
  site_entries result{{gt_conf_entry, gt_conf_percentile_entry}};
  return result;
}

void LevelGenotypedSite::null_model_specific_entries() {
  this->gt_conf = 0.;
  this->gt_conf_percentile = 0.;
}
