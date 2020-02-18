#include "simulate/simulate.hpp"
#include "prg/coverage_graph.hpp"
#include "genotype/infer/json_spec/prg_spec.hpp"
#include "genotype/infer/allele_extracter.hpp"

using namespace gram::simulate;

RandomGenotyper::RandomGenotyper(coverage_Graph const& cov_graph){
    this->cov_graph = &cov_graph;
    child_m = build_child_map(cov_graph.par_map); // Required for site invalidation
    genotyped_records.resize(cov_graph.bubble_map.size()); // Pre-allocate one slot for each bubble in the PRG

    auto json_ptr = std::make_shared<Json_Prg>();
    this->json_prg = json_ptr;

    // Genotype each bubble in the PRG, in most nested to less nested order.
    for (auto const& bubble_pair : cov_graph.bubble_map) {
        auto site_ID = bubble_pair.first->get_site_ID();
        auto site_index = siteID_to_index(site_ID);

        auto extracter = AlleleExtracter(bubble_pair.first, bubble_pair.second, genotyped_records);
    }
}

void gram::commands::simulate::run(SimulateParams const& parameters){

}
