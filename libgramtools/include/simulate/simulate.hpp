#ifndef GRAMTOOLS_SIMULATE_HPP
#define GRAMTOOLS_SIMULATE_HPP

#include "parameters.hpp"
#include "genotype/infer/interfaces.hpp"
#include "common/random.hpp"

using namespace gram::genotype::infer;

namespace gram::simulate {
    class RandomGenotypedSite : public GenotypedSite {
    public:
        RandomGenotypedSite();
        void add_model_specific_JSON(JSON &input_json) override {}
    };

    using Seed = uint32_t;

    class RandomGenotyper : public Genotyper {
    public:
        /**
         * The genotyping process is the same in form to `gram::genotype::infer::LevelGenotyper`
         * except that genotype is randomly assigned among the list of alleles.
         */
        RandomGenotyper(coverage_Graph const &cov_graph);
    };

    gt_site_ptr make_randomly_genotyped_site(RandomGenerator const *const rand, allele_vector const& alleles);
}

namespace gram::commands::simulate {
    void run(SimulateParams const &parameters);
}

#endif //GRAMTOOLS_SIMULATE_HPP
