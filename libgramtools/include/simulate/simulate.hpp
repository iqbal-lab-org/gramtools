#ifndef GRAMTOOLS_SIMULATE_HPP
#define GRAMTOOLS_SIMULATE_HPP

#include "parameters.hpp"
#include "genotype/infer/interfaces.hpp"
using namespace gram::genotype::infer;

namespace gram::simulate {
    class RandomGenotypedSite : public GenotypedSite {
    public:
        void add_model_specific_JSON(JSON &input_json) override {}
    };

    class RandomGenotyper : public Genotyper {
    public:
        RandomGenotyper(coverage_Graph const& cov_graph);
    };
}

namespace gram::commands::simulate {
    void run(SimulateParams const &parameters);
}

#endif //GRAMTOOLS_SIMULATE_HPP
