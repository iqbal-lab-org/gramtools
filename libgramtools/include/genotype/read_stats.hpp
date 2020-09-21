/** @file
 *  Compute read-related statistics.
 *  For inferring a personalised reference, we need estimation of base-level
 * error-rate, and estimation of expected coverage. This module allows for
 * recording that, as well as some other usable metrics, such as max read length
 * and number of sites with no coverage.
 */
#include "genotype/quasimap/coverage/types.hpp"
#include "prg/types.hpp"

#ifndef GRAMTOOLS_READSTATS_HPP
#define GRAMTOOLS_READSTATS_HPP

#define NUM_READS_USED 10000

namespace gram {

class AbstractReadStats {
 public:
  AbstractReadStats()
      : mean_cov_depth(-1),
        variance_cov_depth(-1),
        num_sites_noCov(0),
        num_sites_total(-1) {}
  using allele_and_cov = std::pair<Allele, CovCount>;
  virtual ~AbstractReadStats(){};

  /**
   * Extracts single allele with max. coverage from `start_node` and `end_node`
   * delimiting a variant site in the coverage graph.
   */
  virtual allele_and_cov extract_max_coverage_allele(
      SitesGroupedAlleleCounts const& gped_covs, covG_ptr start_node,
      covG_ptr end_node) = 0;

  /**
   * Compute the depth of coverage using recorded coverage of reads over variant
   * sites after `quasimap`.
   * @param coverage gram::Coverage containing read coverage over variant sites.
   */
  void compute_coverage_depth(Coverage const& coverage,
                              coverage_Graph const& cov_graph);

  double const& get_mean_cov() const { return mean_cov_depth; }
  double const& get_var_cov() const { return variance_cov_depth; }
  std::size_t const& get_num_sites_noCov() const { return num_sites_noCov; }
  std::size_t const& get_num_sites_total() const { return num_sites_total; }

 protected:
  double mean_cov_depth;  // Mean of total coverages at each variant site
  double variance_cov_depth;
  std::size_t num_sites_noCov;
  std::size_t num_sites_total;
};

/**
 * Class for recording read-related statistics.
 */
class ReadStats : public AbstractReadStats {
 public:
  // Default constructor: -1 initialisation to signal that attribute has not
  // been computed.
  ReadStats()
      : no_qual_reads(-1),
        max_read_length(0),
        num_bases_processed(-1),
        mean_pb_error(-1) {}

  /**
   * From a file
   */
  void compute_base_error_rate(const std::string& reads_fpath);

  /**
   * From random access memory
   */
  void compute_base_error_rate(GenomicRead_vector const& reads);

  /**
   * Compute probability of erroneous base from base Phred scores.
   */
  void process_read_perbase_error_rates(AbstractGenomicReadIterator& reads_it);

  using haplogroup_cov = std::pair<AlleleId, CovCount>;

  static haplogroup_cov get_max_cov_haplogroup(
      GroupedAlleleCounts const& gped_cov);

  allele_and_cov extract_max_coverage_allele(
      SitesGroupedAlleleCounts const& gped_covs, covG_ptr start_node,
      covG_ptr end_node) override;

  void serialise(const std::string& json_output_fpath);

  double const& get_mean_pb_error() const { return mean_pb_error; }
  int64_t const& get_num_bases_processed() const { return num_bases_processed; }
  std::size_t const& get_max_read_len() const { return max_read_length; }
  int64_t const& get_num_no_qual_reads() const { return no_qual_reads; }

 private:
  double mean_pb_error;  // Pb sequencing error rate
  int64_t no_qual_reads;
  std::size_t max_read_length;
  int64_t num_bases_processed;
};

}  // namespace gram

#endif  // GRAMTOOLS_READSTATS_HPP
