/** @file
 * Code responsible for inducing input sequences in a genome graph.
 * For each sequence it produces a set of genotyped sites representing its
 * traversal (if one was found- else throws).
 */
#ifndef GMTOOLS_SIMU_INDUCE_GTS_HPP
#define GMTOOLS_SIMU_INDUCE_GTS_HPP

#include <memory>
#include <vector>

#include "genotype/infer/level_genotyping/site.hpp"
#include "genotype/infer/output_specs/fields.hpp"
#include "prg/coverage_graph.hpp"

namespace gram::simulate {
using namespace genotype::infer;
class SimulatedSite : public LevelGenotypedSite {
 public:
  SimulatedSite() = default;
  site_entries get_model_specific_entries() override { return {}; }
  void null_model_specific_entries() override {}
};

class NodeThread;
using nt_ptr = std::shared_ptr<const NodeThread>;
using nt_ptr_v = std::vector<nt_ptr>;

class NoEndpoints : public std::runtime_error {
  using std::runtime_error::runtime_error;
};
class TooManyEndpoints : public std::runtime_error {
  using std::runtime_error::runtime_error;
};

class NodeThread : public std::enable_shared_from_this<NodeThread const> {
 public:
  explicit NodeThread(nt_ptr const& input_parent, covG_ptr input_prg_node,
                      int input_offset)
      : parent(input_parent), prg_node(input_prg_node), offset(input_offset) {}

  ~NodeThread() = default;

  NodeThread(NodeThread const& other) = delete;
  NodeThread(NodeThread const&& other) = delete;
  NodeThread& operator=(NodeThread const& other) = delete;
  NodeThread& operator=(NodeThread const&& other) = delete;

  nt_ptr const& get_parent() const { return parent; }
  covG_ptr const& get_prg_node() const { return prg_node; }
  int get_offset() const { return offset; }
  bool has_next() const { return prg_node->get_num_edges() > 0; }
  void visit(nt_ptr_v& to_visit, std::string const& sequence) const;

 private:
  nt_ptr parent;
  covG_ptr prg_node;
  int offset;
};

/**
 * Makes null genotype calls at all sites in a coverage graph
 */
gt_sites make_nulled_sites(coverage_Graph const& input_prg);

/**
 * Finds all occurrences of `sequence` in the graph
 * @return a vector of endpoints, which are vectors of `NodeThread`s
 */
nt_ptr_v thread_sequence(covG_ptr root, std::string const& sequence);
/**
 * Returns a single endpoint if there is one, else throws.
 * If there are multiple endpoints:
 *  - the first element of the returned pair is set to true.
 *  - the second element is the endpoint that consumed the most of the input
 * sequence
 */
std::pair<bool, nt_ptr> get_single_endpoint(nt_ptr_v const& endpoints,
                                            std::string const& seq_id,
                                            bool const no_ambiguous);
/**
 * Populates traversed @param sites with genotyping information: ref allele,
 * called allele, and sets AMBIG filter if @param has_ambiguity is true.
 */
void apply_genotypes(nt_ptr const& end_point, bool const has_ambiguity,
                     gt_sites const& sites);
gt_sites induce_genotypes_one_seq(gt_sites const& template_sites,
                                  coverage_Graph const& input_prg,
                                  std::string const& sequence,
                                  std::string const& seq_id);
}  // namespace gram::simulate

#endif  // GMTOOLS_SIMU_INDUCE_GTS_HPP
