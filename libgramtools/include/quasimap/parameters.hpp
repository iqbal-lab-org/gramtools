#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/variant/variant.hpp>
#include <boost/variant/get.hpp>
#include <boost/filesystem.hpp>

#include "common/parameters.hpp"


namespace po = boost::program_options;
namespace fs = boost::filesystem;


#ifndef GRAMTOOLS_QUASIMAP_PARAMETERS_HPP
#define GRAMTOOLS_QUASIMAP_PARAMETERS_HPP

namespace gram::commands::quasimap {
    Parameters parse_parameters(po::variables_map &vm,
                                const po::parsed_options &parsed);
}

#endif //GRAMTOOLS_QUASIMAP_PARAMETERS_HPP
