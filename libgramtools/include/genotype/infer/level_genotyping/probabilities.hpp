#ifndef GRAMTOOLS_PROBS_HPP
#define GRAMTOOLS_PROBS_HPP

#include <map>
#include <memory>
#include <vector>

#include "genotype/quasimap/coverage/types.hpp"

namespace gram::genotype::infer::probabilities {
using params = std::vector<double>;
using memoised_params = std::map<params, double>;

class AbstractPmf {
 protected:
  AbstractPmf() = default;
  memoised_params probs;  // Memoised probabilities
  virtual double compute_prob(params const& query) const = 0;

 public:
  virtual ~AbstractPmf() = default;
  double operator()(params const& query);
  memoised_params const& get_probs() const { return probs; }
};

class PoissonLogPmf : public AbstractPmf {
  double lambda;
  double compute_prob(params const& query) const override;

 public:
  PoissonLogPmf() : lambda(0) {}
  explicit PoissonLogPmf(params const& parameterisation);
};

class NegBinomLogPmf : public AbstractPmf {
  // k: number of successes; p: probability of success
  double k, p;
  double compute_prob(params const& query) const override;

 public:
  explicit NegBinomLogPmf(params const& parameterisation);
};

using pmf_ptr = std::shared_ptr<AbstractPmf>;
struct DataParams {
  double mean_cov{-1};
  double mean_pb_error{-1};
  double num_successes{-1};
  double success_prob{-1};
  DataParams() = default;
  DataParams(double const input_mean_cov, double const input_mean_pb_error)
      : mean_cov(input_mean_cov), mean_pb_error(input_mean_pb_error) {}
  bool operator==(DataParams const& other) const {
    return mean_cov == other.mean_cov && mean_pb_error == other.mean_pb_error &&
           num_successes == other.num_successes &&
           success_prob == other.success_prob;
  }
};

struct likelihood_related_stats {
  DataParams data_params;
  double log_mean_pb_error, log_zero, log_zero_half_depth, log_no_zero,
      log_no_zero_half_depth;
  CovCount credible_cov_t; /**< minimum coverage count to qualify as actual
                              coverage (per-base)*/
  pmf_ptr pmf_full_depth;
  pmf_ptr pmf_half_depth;
};

std::ostream& operator<<(std::ostream& os,
                         const likelihood_related_stats& l_stats);

}  // namespace gram::genotype::infer::probabilities

#endif  // GRAMTOOLS_PROBS_HPP
