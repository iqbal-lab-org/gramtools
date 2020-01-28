#ifndef TEST_SRC_COMMON
#define TEST_SRC_COMMON

#include "src_common/generate_prg.hpp"
#include "kmer_index/build.hpp"
#include "genotype/quasimap/coverage/common.hpp"
#include "genotype/read_stats.hpp"
#include "prg/coverage_graph.hpp"

using namespace gram;

using prg_positions = std::vector<std::size_t>;
using covG_ptrPair = std::pair<covG_ptr, covG_ptr>;

/**
 * Given a `cov_graph` and a set of positions in the PRG string,
 * returns the coverage of each node in the coverage graph corresponding to each position.
 *
 * Useful for testing per base coverage recordings.
 */
gram::SitePbCoverage collect_coverage(coverage_Graph const& cov_graph, prg_positions positions);

/**
 * Given a map of all bubbles and a `siteID` of interest, returns the pair of `covG_ptr` corresponding
 * to the start and end nodes of the site.
 */
covG_ptrPair get_bubble_nodes(covG_ptr_map bubble_map, Marker site_ID);

/**
 * Builds a coverage graph, fm-index and kmer index from a PRG string.
 * Particularly useful in `genotype` steps: quasimap and infer.
 */
class prg_setup{
public:
    PRG_Info prg_info;
    Coverage coverage;
    Parameters parameters;
    KmerIndex kmer_index;
    ReadStats read_stats;

    explicit prg_setup() {};

    /**
     * Sets up a 'legacy'-style PRG string, with no nesting
     */
    void setup_numbered_prg(std::string raw_prg,
                            Sequences kmers){
        auto encoded_prg = encode_prg(raw_prg);
        internal_setup(encoded_prg, kmers);
    }

    /**
     * The bracketed format allows unambiguously encoding nested PRG strings.
     */
    void setup_bracketed_prg(std::string raw_prg,
                             Sequences kmers){
        auto encoded_prg = prg_string_to_ints(raw_prg);
        internal_setup(encoded_prg, kmers);
    }

    /**
     * Maps reads and populates ReadStats instance from the raw reads and the mapped instances
     */
    void quasimap_reads( GenomicRead_vector const& reads );

private:
    void internal_setup(marker_vec encoded_prg, Sequences kmers);
};


#endif //TEST_SRC_COMMON

