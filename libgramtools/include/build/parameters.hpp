/**
 * @file
 * Command-line argument processing for `build` command.
 */

#ifndef GRAMTOOLS_BUILD_PARAMETERS_HPP
#define GRAMTOOLS_BUILD_PARAMETERS_HPP

#include "common/parameters.hpp"

namespace gram {
class BuildParams : public CommonParameters {
 public:
  std::string sdsl_memory_log_fpath;
  uint32_t max_read_size;
  bool all_kmers_flag;
  std::string fasta_ref;
};

namespace commands::build {
/**
 * Parse command-line parameters to `build` option.
 * Options include kmer size for index and whether to build all kmers or only
 * those needed from the prg.
 * @return `Parameters`
 */
BuildParams parse_parameters(po::variables_map &vm,
                             const po::parsed_options &parsed);
}  // namespace commands::build
}  // namespace gram

#endif  // GRAMTOOLS_BUILD_PARAMETERS_HPP
