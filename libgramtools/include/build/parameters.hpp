/**
 * @file
 * Command-line argument processing for `build` command.
 */
namespace po = boost::program_options;


#ifndef GRAMTOOLS_BUILD_PARAMETERS_HPP
#define GRAMTOOLS_BUILD_PARAMETERS_HPP

namespace gram::commands::build {
    /**
     * Parse command-line parameters to `build` option.
     * Options include kmer size for index and whether to build all kmers or only those needed from the prg.
     * @return `Parameters`
     */
    Parameters parse_parameters(po::variables_map &vm, const po::parsed_options &parsed);
}

#endif //GRAMTOOLS_BUILD_PARAMETERS_HPP
