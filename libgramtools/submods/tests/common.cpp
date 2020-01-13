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

covG_ptrPair get_bubble_nodes(covG_ptr_map bubble_map, Marker site_ID){
    ensure_is_site_marker(site_ID);
    for (auto const& bubble : bubble_map){
        if (bubble.first->get_site_ID() == site_ID) return bubble;
    }
    throw std::invalid_argument("The provided site ID was not found in the map of PRG bubbles.");
}
