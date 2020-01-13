#ifndef TEST_SRC_COMMON
#define TEST_SRC_COMMON

#include "genotype/quasimap/coverage/types.hpp"
#include "prg/coverage_graph.hpp"

using prg_positions = std::vector<std::size_t>;
using covG_ptrPair = std::pair<covG_ptr, covG_ptr>;

/**
 * Given a `cov_graph` and a set of positions in the PRG string,
 * returns the coverage of each node in the coverage graph corresponding to each position.
 *
 * Useful for testing per base coverage recordings.
 */
gram::AlleleCoverage collect_coverage(coverage_Graph const& cov_graph, prg_positions positions);

/**
 * Given a map of all bubbles and a `siteID` of interest, returns the pair of `covG_ptr` corresponding
 * to the start and end nodes of the site.
 */
covG_ptrPair get_bubble_nodes(covG_ptr_map bubble_map, Marker site_ID);

#endif //TEST_SRC_COMMON

