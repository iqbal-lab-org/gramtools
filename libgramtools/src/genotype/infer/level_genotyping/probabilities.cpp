#include "genotype/infer/level_genotyping/probabilities.hpp"

#include <assert.h>

#include <cmath>

namespace gram::genotype::infer::probabilities {
double AbstractPmf::operator()(params const& query) {
  if (probs.find(query) != probs.end())
    return probs.at(query);
  else {
    auto prob = compute_prob(query);
    probs.insert(std::pair<params, double>(query, prob));
    return prob;
  }
}

double PoissonLogPmf::compute_prob(params const& query) const {
  assert(query.size() == 1);
  auto cov{query.at(0)};
  return (-1 * lambda + cov * log(lambda) - lgamma(cov + 1));
}

PoissonLogPmf::PoissonLogPmf(params const& parameterisation)
    : lambda(parameterisation[0]) {
  operator()(params{0});
}

NegBinomLogPmf::NegBinomLogPmf(params const& parameterisation)
    : k(parameterisation[0]), p(parameterisation[1]), AbstractPmf() {
  operator()(params{0});
}

double NegBinomLogPmf::compute_prob(params const& query) const {
  assert(query.size() == 1);
  auto cov{query.at(0)};
  return (lgamma(k + cov) - lgamma(cov + 1) - lgamma(k) + k * log(p) +
          cov * log(1 - p));
}

std::ostream& operator<<(std::ostream& os,
                         const likelihood_related_stats& l_stats) {
  char buf[1024];
  std::sprintf(buf,
               "Model params: \nmean cov: %f\n"
               "mean per-base error: %f\nnum successes: %f\n"
               "prob of success: %f \n"
               "log_prob_zero_cov: %f \n"
               "log_prob_nonzero_cov: %f\n",
               l_stats.data_params.mean_cov, l_stats.data_params.mean_pb_error,
               l_stats.data_params.num_successes,
               l_stats.data_params.success_prob, l_stats.log_zero,
               l_stats.log_no_zero);
  os << buf;
  return os;
}
}  // namespace gram::genotype::infer::probabilities
