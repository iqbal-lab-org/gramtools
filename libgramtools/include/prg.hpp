#include <cctype>
#include <cstdlib>
#include <vector>
#include <string>
#include <tuple>

#include "masks.hpp"
#include "fm_index.hpp"
#include "ranks.hpp"


#ifndef GRAMTOOLS_PRG_HPP
#define GRAMTOOLS_PRG_HPP

struct PRG_Info {
    FM_Index fm_index;
    DNA_Rank dna_rank;
    SitesMask sites_mask;
    AlleleMask allele_mask;
    uint64_t max_alphabet_num;
};

uint64_t max_alphabet_num(const std::string &prg_raw);

std::vector<int> generate_allele_mask(const std::string &prg_raw);

#endif //GRAMTOOLS_PRG_HPP
