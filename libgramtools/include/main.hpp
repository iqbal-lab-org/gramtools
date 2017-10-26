#include "parameters.hpp"


#ifndef GRAMTOOLS_MAIN_HPP
#define GRAMTOOLS_MAIN_HPP

// Parameters parse_command_line_parameters(int argc, const char *const *argv);
std::pair<Parameters, Commands> parse_command_line_parameters(int argc, const char *const *argv);

void build(const Parameters &parameters);

void quasimap(const Parameters &params);

PRG_Info load_prg_info(const Parameters &params);

#endif //GRAMTOOLS_MAIN_HPP
