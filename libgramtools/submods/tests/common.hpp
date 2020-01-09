#ifndef TEST_SRC_COMMON
#define TEST_SRC_COMMON

#include "genotype/quasimap/coverage/types.hpp"
#include "prg/coverage_graph.hpp"

using prg_positions = std::vector<std::size_t>;

/**
 * Given a @param cov_graph and a set of positions in the PRG string,
 * returns the coverage of each node in the coverage graph corresponding to each position.
 *
 * Useful for testing per base coverage recordings.
 */
gram::AlleleCoverage collect_coverage(coverage_Graph const& cov_graph, prg_positions positions);

#endif //TEST_SRC_COMMON

