#ifndef LVLGT_MODEL
#define LVLGT_MODEL

#include "genotype/parameters.hpp"
#include "probabilities.hpp"
#include "site.hpp"

using namespace gram;

namespace gram::genotype::infer {
class UnsupportedPloidy : public std::runtime_error {
  using std::runtime_error::runtime_error;
};

class IncorrectGenotyping : public std::runtime_error {
  using std::runtime_error::runtime_error;
};

using multiplicities = std::vector<bool>;
using likelihood_map = std::multimap<double, GtypedIndices, std::greater<>>;
using memoised_coverages = std::map<AlleleIds, allele_coverages>;
using CovPair = std::pair<double, double>;

using namespace probabilities;

struct ModelData {
  allele_vector const input_alleles;
  GroupedAlleleCounts const gp_counts;
  Ploidy ploidy;
  likelihood_related_stats const *l_stats;
  bool debug = false;

  ModelData() : gp_counts() {}

  ModelData(allele_vector const &input_alleles,
            GroupedAlleleCounts const &gp_counts, Ploidy ploidy,
            likelihood_related_stats const *l_stats, bool debug = false)
      : input_alleles(input_alleles),
        gp_counts(gp_counts),
        ploidy(ploidy),
        l_stats(l_stats),
        debug(debug){};
};

/**
Genotyping model using:
  * coverage equivalence-classes
  * alternative alleles all at same nesting level
  * genotype confidence using likelihood ratios
  * invalidation of nested bubbles
*/
class LevelGenotyperModel : GenotypingModel {
  using site_ptr = std::shared_ptr<LevelGenotypedSite>;
  ModelData data;

  // Computed at construction time
  PerAlleleCoverage haploid_allele_coverages;   /**< Coverage counts compatible
                                                   with single alleles */
  PerAlleleCoverage singleton_allele_coverages; /**< Coverage counts unique to
                                                   single alleles */
  memoised_coverages computed_coverages;
  std::size_t total_coverage;

  // Computed at run time
  likelihood_map likelihoods;  // Stores highest likelihoods first
  site_ptr genotyped_site;     // What the class will build

 public:
  LevelGenotyperModel() = default;
  explicit LevelGenotyperModel(ModelData &input_data);

  bool ignore_ref_allele() const {
    return !data.input_alleles.at(0).nesting_consistent;
  }

  // Constructor for testing
  LevelGenotyperModel(likelihood_related_stats const &input_l_stats,
                      PerAlleleCoverage const &input_covs,
                      likelihood_map const &input_likelihoods);

  /*_______Preparations______*/
  std::size_t count_total_coverage(GroupedAlleleCounts const &gp_counts);

  std::vector<bool> get_haplogroup_multiplicities(
      allele_vector const &input_alleles);

  void set_haploid_coverages(GroupedAlleleCounts const &input_gp_counts,
                             AlleleId num_haplogroups);
  /**
   * Alleles with no sequence correspond to direct deletions.
   * In this case they get assigned coverage by this function,
   * using the grouped allele coverages, as if they had a single base.
   */
  void assign_coverage_to_empty_alleles(allele_vector &input_alleles);

  /*_______Likelihoods______*/
  /**
   * Counts the number of positions in an allele with coverage below threshold
   * `credible_cov_t`. This threshold is the coverage at which true coverage is
   * more likely than erroneous (sequencing error-based) coverage.
   */
  double fraction_noncredible_positions(Allele const &allele);

  /**
   * Computes log-likelihood of allelic coverage and stores it.
   * Handles haploid and diploid flexibly.
   */
  void add_likelihood(allele_vector const &alleles,
                      double const &incompatible_coverage,
                      GtypedIndices const &allele_indices);

  /**
   * Haploid genotype likelihood
   */
  void compute_haploid_log_likelihoods(allele_vector const &input_alleles);

  /**
   * Diploid homozygous
   */
  void compute_homozygous_log_likelihoods(
      allele_vector const &input_alleles,
      multiplicities const &haplogroup_multiplicities);

  /**
   * Diploid. Because of the large possible number of diploid combinations,
   * (eg for 10 alleles, 45), we only consider for combination those alleles
   * that have at least one unit of coverage unique to them.
   */
  void compute_heterozygous_log_likelihoods(
      allele_vector const &input_alleles,
      multiplicities const &haplogroup_multiplicities);

  /** For producing the diploid combinations. */
  std::vector<GtypedIndices> get_permutations(GtypedIndices const &indices,
                                              std::size_t subset_size);

  /*_______Coverages______*/
  /**
   * Note: Due to nesting, the alleles can be from the same haplogroup; in which
   * case, they have the same haploid coverage, and they get assigned half of it
   * each.
   */
  std::pair<double, double> compute_diploid_coverage(
      GroupedAlleleCounts const &gp_counts, AlleleIds haplogroups,
      multiplicities const &haplogroup_multiplicities);

  std::pair<double, double> diploid_cov_same_haplogroup(
      AlleleIds const &haplogroups);
  std::pair<double, double> diploid_cov_different_haplogroup(
      GroupedAlleleCounts const &gp_counts, AlleleIds const &ids,
      multiplicities const &hap_mults);

  /*_______Make result______*/
  void CallGenotype(allele_vector const &input_alleles,
                    multiplicities hap_mults, Ploidy const ploidy);

  /**
   * If coverage differences between chosen allele(s) and next best allele(s)
   * is small, adds the next best alleles for consideration by parent sites,
   * if any.
   * This propagates uncertainty upwards, reducing overconfidence.
   */
  void add_next_best_alleles(allele_vector const &input_alleles,
                             GtypedIndices const &chosen_gt,
                             GtypedIndices const &next_best_gt);

  /**
   * If there is no coverage difference between chosen allele(s) and next best
   * allele(s), adds all these alleles for consideration by parent sites,
   * if any.
   * This propagates uncertainty upwards, reducing overconfidence.
   */
  void add_all_best_alleles(allele_vector const &input_alleles,
                            GtypedIndices const &chosen_gt,
                            GtypedIndices const &next_best_gt);

  /**
   * Picks best likelihood conditional on the genotype not being
   * nesting inconsistent, to ensure model cannot choose a parent site allele
   * inconsistent with child site calls.
   */
  static likelihood_map::const_iterator ChooseMaxLikelihood(
      likelihood_map const &likelihoods, allele_vector const &alleles);

  AlleleIds get_haplogroups(allele_vector const &alleles,
                            GtypedIndices const &gtype) const;

  /**
   * Express genotypes as relative to chosen alleles.
   * For eg, {0, 2, 4} in the original set of possible alleles goes to {0, 1, 2}
   * in the 3 called alleles.
   */
  GtypedIndices rescale_genotypes(GtypedIndices const &genotypes);

  // Getters
  PerAlleleCoverage const &get_haploid_covs() const {
    return haploid_allele_coverages;
  }
  PerAlleleCoverage const &get_singleton_covs() const {
    return singleton_allele_coverages;
  }
  likelihood_map const &get_likelihoods() const { return likelihoods; }

  gt_site_ptr get_site() override {
    return std::static_pointer_cast<gt_site>(genotyped_site);
  }
  gtype_information get_site_gtype_info() {
    return genotyped_site->get_all_gtype_info();
  }
  // For use by GCP library
  double get_genotype_confidence() { return genotyped_site->get_gt_conf(); }
};
}  // namespace gram::genotype::infer

#endif  // LVLGT_MODEL
