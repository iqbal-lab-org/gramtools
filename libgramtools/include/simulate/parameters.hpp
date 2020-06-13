#ifndef SIMU_PARAMETERS_HPP
#define SIMU_PARAMETERS_HPP

#include "common/parameters.hpp"

namespace gram {

class SimulateParams : public CommonParameters {
 public:
  std::string json_out_fpath, fasta_out_fpath;
  std::string sample_id;
  uint64_t max_num_paths;
  std::string input_sequences_fpath;
};

namespace commands::simulate {
SimulateParams parse_parameters(po::variables_map &vm,
                                const po::parsed_options &parsed);
}
}  // namespace gram

#endif  // SIMU_PARAMETERS_HPP
