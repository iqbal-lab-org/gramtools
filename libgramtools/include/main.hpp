/**
 * @file
 * Entry point for gramtools backend (`gram`).
 */
#include "common/parameters.hpp"


#ifndef GRAMTOOLS_MAIN_HPP
#define GRAMTOOLS_MAIN_HPP

namespace gram {
    std::pair<Parameters, Commands> parse_command_line_parameters(int argc, const char *const *argv);
}

#endif //GRAMTOOLS_MAIN_HPP
