namespace po = boost::program_options;


#ifndef GRAMTOOLS_BUILD_PARAMETERS_HPP
#define GRAMTOOLS_BUILD_PARAMETERS_HPP

namespace gram::commands::build {
    Parameters parse_parameters(po::variables_map &vm, const po::parsed_options &parsed);
}

#endif //GRAMTOOLS_BUILD_PARAMETERS_HPP
