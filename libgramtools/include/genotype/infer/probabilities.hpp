#ifndef GRAMTOOLS_PROBS_HPP
#define GRAMTOOLS_PROBS_HPP

#include <vector>
#include <map>

namespace gram::genotype::infer::probabilities{
    using params = std::vector<float>;
    using memoised_params = std::map<params, float>;

    class AbstractPmf{
    protected:
        AbstractPmf() = default;
        AbstractPmf(params const& parameterisation) : parameterisation(parameterisation) {}
        params const parameterisation; // Fixed values of the pmf parameters
        memoised_params probs; // Memoised probabilities
        virtual float compute_prob(params const& query) const = 0;

    public:
        virtual ~AbstractPmf(){}
        float operator()(params const& query);
        memoised_params const& get_probs() const { return probs; }
    };

    class PoissonLogPmf : public AbstractPmf{
        float const lambda;
        float compute_prob(params const& query) const override;
    public:
        PoissonLogPmf() : lambda(0) {}
        PoissonLogPmf(params const& parameterisation);
    };
}

#endif //GRAMTOOLS_PROBS_HPP
