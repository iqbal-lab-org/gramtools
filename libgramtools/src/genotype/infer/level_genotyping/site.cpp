#include "genotype/infer/level_genotyping/site.hpp"
#include "genotype/infer/output_specs/fields.hpp"

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

header_vec LevelGenotypedSite::site_model_specific_entries() {
  return header_vec{
      vcf_meta_info_line{
          "FORMAT", "GT_CONF",
          "Genotype confidence as "
          "likelihood ratio of called and next most likely genotype.",
          "1", "Float"},
      vcf_meta_info_line{"FORMAT", "GT_CONF_PERCENTILE",
                         "Percent of calls expected to have lower GT_CONF", "1",
                         "Float"}};
}

site_entries LevelGenotypedSite::get_model_specific_entries() {
  auto specs = site_model_specific_entries();
  site_entry<double> gc_entry(specs.at(0));
  gc_entry.vals.emplace_back(gt_conf);

  site_entry<double> gcp_entry(specs.at(1));
  gcp_entry.vals.emplace_back(gt_conf_percentile);

  return site_entries{{gc_entry, gcp_entry}};
}

void LevelGenotypedSite::null_model_specific_entries() {
  this->gt_conf = 0.;
  this->gt_conf_percentile = 0.;
}
