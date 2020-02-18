#include "simulate/simulate.hpp"
#include "prg/coverage_graph.hpp"
#include "genotype/infer/json_spec/prg_spec.hpp"
#include "genotype/infer/allele_extracter.hpp"
#include "genotype/infer/personalised_reference.hpp"

using namespace gram::genotype;
using namespace gram::simulate;

namespace gram::simulate {

    RandomGenotyper::RandomGenotyper(coverage_Graph const &cov_graph, Seed const &seed) {
        this->rand = RandomInclusiveInt(seed);
        this->cov_graph = &cov_graph;
        child_m = build_child_map(cov_graph.par_map); // Required for site invalidation
        genotyped_records.resize(cov_graph.bubble_map.size()); // Pre-allocate one slot for each bubble in the PRG

        auto json_ptr = std::make_shared<Simulated_Json>();
        this->json_prg = json_ptr;

        // Genotype each bubble in the PRG, in most nested to less nested order.
        for (auto const &bubble_pair : cov_graph.bubble_map) {
            auto site_ID = bubble_pair.first->get_site_ID();
            auto site_index = siteID_to_index(site_ID);

            auto extracter = AlleleExtracter(bubble_pair.first, bubble_pair.second, genotyped_records);
            auto genotyped_site = make_randomly_genotyped_site(&rand, extracter.get_alleles());

            genotyped_records.at(site_index) = genotyped_site;
            // Line below is so that when allele extraction occurs and jumps through a previously
            // genotyped site, it knows where in the graph to resume from.
            genotyped_records.at(site_index)->set_site_end_node(bubble_pair.second);

            run_invalidation_process(genotyped_site, site_ID);
        }
    }

    gt_site_ptr make_randomly_genotyped_site(RandomGenerator const *const rand, allele_vector const &alleles) {
        allele_vector picked_alleles{alleles.begin(), alleles.begin() + 1};  // Always pick REF
        auto picked_index = rand->generate(0, alleles.size() - 1);

        if (picked_index != 0) {
            picked_alleles.push_back(alleles.at(picked_index));
            picked_index = 1;
        }

        auto result = std::make_shared<RandomGenotypedSite>();
        result->set_genotype(GtypedIndices{picked_index});
        result->set_alleles(picked_alleles);
        return result;
    }
}

void gram::commands::simulate::run(SimulateParams const& parameters){

}
