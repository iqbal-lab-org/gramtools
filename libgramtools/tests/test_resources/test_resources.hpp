#ifndef TEST_SRC_COMMON
#define TEST_SRC_COMMON

#include "build/kmer_index/build.hpp"
#include "genotype/parameters.hpp"
#include "genotype/quasimap/quasimap.hpp"
#include "genotype/read_stats.hpp"
#include "types.hpp"

using namespace gram;

/**
 * Given a `cov_graph` and a set of positions in the PRG string,
 * returns the coverage of each node in the coverage graph corresponding to each
 * position.
 *
 * Useful for testing per base coverage recordings.
 */
gram::SitePbCoverage collect_coverage(coverage_Graph const& cov_graph,
                                      prg_positions positions);

/**
 * Builds a coverage graph, fm-index and kmer index from a PRG string.
 * Particularly useful in `genotype` steps: quasimap and infer.
 */
class prg_setup {
 public:
  PRG_Info prg_info;
  Coverage coverage;
  GenotypeParams parameters;
  KmerIndex kmer_index;
  ReadStats read_stats;
  gram::QuasimapReadsStats quasimap_stats;

  explicit prg_setup(){};

  /**
   * Sets up a 'legacy'-style PRG string, with no nesting
   */
  void setup_numbered_prg(std::string raw_prg, uint32_t kmer_size = 2) {
    auto kmers_set = gram::generate_all_kmers(kmer_size);
    Sequences kmers{kmers_set.begin(), kmers_set.end()};
    auto encoded_prg = encode_prg(raw_prg);
    internal_setup(encoded_prg, kmers);
  }

  /**
   * The bracketed format allows unambiguously encoding nested PRG strings.
   */
  void setup_bracketed_prg(std::string raw_prg, uint32_t kmer_size = 2) {
    auto kmers_set = gram::generate_all_kmers(kmer_size);
    Sequences kmers{kmers_set.begin(), kmers_set.end()};
    auto encoded_prg = prg_string_to_ints(raw_prg);
    internal_setup(encoded_prg, kmers);
  }

  /**
   * Maps reads and populates ReadStats instance from the raw reads and the
   * mapped instances
   */
  void quasimap_reads(GenomicRead_vector const& reads);

 private:
  void internal_setup(marker_vec encoded_prg, Sequences kmers);
};

#endif  // TEST_SRC_COMMON
