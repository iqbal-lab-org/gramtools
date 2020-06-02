#ifndef GRAMTOOLS_PROBS_HPP
#define GRAMTOOLS_PROBS_HPP

#include <vector>
#include <map>
#include <memory>

#include "genotype/quasimap/coverage/types.hpp"

namespace gram::genotype::infer::probabilities{
    using params = std::vector<double>;
    using memoised_params = std::map<params, double>;


    class AbstractPmf{
    protected:
        AbstractPmf() = default;
        memoised_params probs; // Memoised probabilities
        virtual double compute_prob(params const& query) const = 0;

    public:
        virtual ~AbstractPmf() = default;
        double operator()(params const& query);
        memoised_params const& get_probs() const { return probs; }
    };


    class PoissonLogPmf : public AbstractPmf{
        double lambda;
        double compute_prob(params const& query) const override;
    public:
        PoissonLogPmf() : lambda(0) {}
        explicit PoissonLogPmf(params const& parameterisation);
    };

    class NegBinomLogPmf : public AbstractPmf{
        // k: number of successes; p: probability of success
        double k, p;
        double compute_prob(params const& query) const override;
    public:
        explicit NegBinomLogPmf(params const& parameterisation);
    };


    using pmf_ptr = std::shared_ptr<AbstractPmf>;
    struct likelihood_related_stats {
        params data_params;
        double log_mean_pb_error,
                log_zero, log_zero_half_depth,
                log_no_zero, log_no_zero_half_depth;
        CovCount credible_cov_t; /**< minimum coverage count to qualify as actual coverage (per-base)*/
        pmf_ptr pmf_full_depth;
        pmf_ptr pmf_half_depth;
    };
}

#endif //GRAMTOOLS_PROBS_HPP
