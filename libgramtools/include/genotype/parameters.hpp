/**
 * @file
 * Command-line argument processing for `quasimap` command.
 */
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/variant/variant.hpp>
#include <boost/variant/get.hpp>
#include <boost/filesystem.hpp>

#include "common/parameters.hpp"

namespace po = boost::program_options;

#ifndef GRAMTOOLS_QUASIMAP_PARAMETERS_HPP
#define GRAMTOOLS_QUASIMAP_PARAMETERS_HPP

namespace gram::commands::genotype {
    /**
     * Parse command line parameters.
     * A directory containing the information necessary for vBWT mapping to the prg must be passed.
     * This contains files produced in the `build` stage. They are then loaded during `quasimap`ping.
     * @see gram::commands::build::run()
     */
    Parameters parse_parameters(po::variables_map &vm,
                                const po::parsed_options &parsed);
}

#endif //GRAMTOOLS_QUASIMAP_PARAMETERS_HPP
