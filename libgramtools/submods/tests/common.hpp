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
    void setup(std::string raw_prg,
               Sequences kmers){
        auto encoded_prg = encode_prg(raw_prg);
        internal_setup(encoded_prg, kmers);
    }

    void setup_nested(std::string raw_prg,
                      Sequences kmers){
        auto encoded_prg = prg_string_to_ints(raw_prg);
        internal_setup(encoded_prg, kmers);
    }

    // TODO: calls gram::quasimap_read and populates ReadStats instance
    void quasimap_reads() {;}

private:
    void internal_setup(marker_vec encoded_prg, Sequences kmers){
        size_t kmer_size = kmers.front().size();
        for (auto const& kmer : kmers) assert(kmer_size == kmer.size());

        // TODO: the calls to rank_support setup in `generate_prg_info` get somehow lost when leaving its scope
        // and we need to call `init_support`, or rank_support again, in this scope for it to work
        prg_info = generate_prg_info(encoded_prg);

        sdsl::util::init_support(prg_info.rank_bwt_a, &prg_info.dna_bwt_masks.mask_a);
        sdsl::util::init_support(prg_info.rank_bwt_c, &prg_info.dna_bwt_masks.mask_c);
        sdsl::util::init_support(prg_info.rank_bwt_g, &prg_info.dna_bwt_masks.mask_g);
        sdsl::util::init_support(prg_info.rank_bwt_t, &prg_info.dna_bwt_masks.mask_t);

        sdsl::util::init_support(prg_info.prg_markers_rank, &prg_info.prg_markers_mask);
        sdsl::util::init_support(prg_info.prg_markers_select, &prg_info.prg_markers_mask);

        coverage = coverage::generate::empty_structure(prg_info);

        parameters.kmers_size = kmer_size;
        kmer_index = index_kmers(kmers, parameters.kmers_size, prg_info);
    }
};


#endif //TEST_SRC_COMMON

