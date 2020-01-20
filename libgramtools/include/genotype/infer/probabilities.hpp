#ifndef GRAMTOOLS_PROBS_HPP
#define GRAMTOOLS_PROBS_HPP

#include <vector>
#include <map>

namespace gram::genotype::infer::probabilities{
    using params = std::vector<double>;
    using memoised_params = std::map<params, double>;

    class AbstractPmf{
    protected:
        AbstractPmf() = default;
        memoised_params probs; // Memoised probabilities
        virtual double compute_prob(params const& query) const = 0;

    public:
        virtual ~AbstractPmf(){}
        double operator()(params const& query);
        memoised_params const& get_probs() const { return probs; }
    };

    class PoissonLogPmf : public AbstractPmf{
        double const lambda;
        double compute_prob(params const& query) const override;
    public:
        PoissonLogPmf() : lambda(0) {}
        PoissonLogPmf& operator=(const PoissonLogPmf& other) {}
        PoissonLogPmf(params const& parameterisation);
    };
}

#endif //GRAMTOOLS_PROBS_HPP
