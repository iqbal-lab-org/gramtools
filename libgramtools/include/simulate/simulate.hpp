#ifndef GRAMTOOLS_SIMULATE_HPP
#define GRAMTOOLS_SIMULATE_HPP

#include "common/random.hpp"
#include "genotype/infer/level_genotyping/runner.hpp"
#include "parameters.hpp"

using namespace gram::genotype::infer;

namespace gram::simulate {
using Seed = uint32_t;

class SimulationGenotyper : public LevelGenotyper {
 public:
  /**
   * The genotyping process is the same in form to
   * `gram::genotype::infer::LevelGenotyper` except that genotype is randomly
   * assigned among the list of alleles.
   */
  SimulationGenotyper(coverage_Graph const &cov_graph);

  /**
   * For taking in directly genotyped sites
   */
  SimulationGenotyper(coverage_Graph const &cov_graph,
                      gt_sites const &input_sites);

  header_vec get_model_specific_headers();
};

/**
 * @param use_ref_allele : if the ref allele is not consistent with the inner
 bubbles, it is still produced by the `AlleleExtracter`. In that case this param
 should be false.
 */
lvlgt_site_ptr make_randomly_genotyped_site(RandomGenerator *const rand,
                                            allele_vector const &alleles);
}  // namespace gram::simulate

namespace gram::commands::simulate {
void run(SimulateParams const &parameters);
}

#endif  // GRAMTOOLS_SIMULATE_HPP
