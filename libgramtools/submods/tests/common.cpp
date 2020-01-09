#include "tests/common.hpp"
using namespace gram;

AlleleCoverage collect_coverage(coverage_Graph const& cov_graph, prg_positions positions){
    AlleleCoverage result(positions.size());
    covG_ptr accessed_node;
    std::size_t index{0};

    for (auto& pos : positions) {
        accessed_node = cov_graph.random_access[pos].node;
        result[index] = accessed_node->get_coverage();
        index++;
    }
    return result;
}

