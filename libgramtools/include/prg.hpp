#include <cctype>
#include <cstdlib>
#include <vector>
#include <string>
#include <tuple>


#ifndef GRAMTOOLS_PRG_HPP
#define GRAMTOOLS_PRG_HPP

uint64_t max_alphabet_num(const std::string &prg_raw);

std::vector<int> generate_allele_mask(const std::string &prg_raw);

#endif //GRAMTOOLS_PRG_HPP
