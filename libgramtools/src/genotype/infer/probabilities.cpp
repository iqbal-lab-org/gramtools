#include <cmath>
#include <assert.h>
#include "genotype/infer/probabilities.hpp"

namespace gram::genotype::infer::probabilities{
    float AbstractPmf::operator()(params const& query){
        if (probs.find(query) != probs.end()) return probs.at(query);
        else {
            auto prob = compute_prob(query);
            probs.insert(std::pair<params, float>(query, prob));
            return prob;
        }
    }


    float PoissonLogPmf::compute_prob(params const& query) const {
        assert(query.size() == 1);
        auto coverage_count{query.at(0)};
        return (-1 * lambda + coverage_count * log(lambda) - lgamma(coverage_count + 1));
    }

    PoissonLogPmf::PoissonLogPmf(params const& parameterisation) :
            AbstractPmf(parameterisation), lambda(parameterisation[0]){
        operator()(params{0}); // Adds Poisson(0) in because is always used in genotyping model
    }
}
