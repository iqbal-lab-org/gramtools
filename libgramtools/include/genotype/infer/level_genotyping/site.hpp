#ifndef LVLGT_SITE
#define LVLGT_SITE

#include <variant>

#include "genotype/infer/interfaces.hpp"

namespace gram::genotype::output_spec {
struct vcf_meta_info_line;
using header_vec = std::vector<vcf_meta_info_line>;
}  // namespace gram::genotype::output_spec

namespace gram::genotype::infer {

class LevelGenotypedSite : public GenotypedSite {
  double gt_conf = 0.; /**< Difference in log likelihood between most likely and
                          next most likely genotype **/
  double gt_conf_percentile =
      0.; /**< Percent of gt_confs that are below this value */
 public:
  ~LevelGenotypedSite() override = default;

  void set_gt_conf(double const& val) { gt_conf = val; }
  double get_gt_conf() const { return gt_conf; }
  void set_gt_conf_percentile(double const& val) { gt_conf_percentile = val; }

  void set_num_haplogroups(std::size_t const& num_haps) {
    num_haplogroups = num_haps;
  }

  /**
   * Produce the haplogroups that have not been genotyped, for use in nested
   * site invalidation.
   */
  AlleleIds const get_nonGenotyped_haplogroups() const;

  AlleleIds const get_all_haplogroups() const {
    assert(num_haplogroups > 0);
    AlleleIds result;
    for (std::size_t idx{0}; idx < num_haplogroups; idx++)
      result.push_back(idx);
    return result;
  }

  static header_vec site_model_specific_entries();
  site_entries get_model_specific_entries() override;
  void null_model_specific_entries() override;
};
}  // namespace gram::genotype::infer

#endif  // LVLGT_SITE
