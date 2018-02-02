#include <fstream>

#include "search/search.hpp"
#include "quasimap/utils.hpp"


uint64_t get_number_of_variant_sites(const PRG_Info &prg_info) {
    auto min_boundary_marker = 5;
    uint64_t numer_of_variant_sites;
    if (prg_info.max_alphabet_num <= 4) {
        numer_of_variant_sites = 0;
    } else {
        numer_of_variant_sites = (prg_info.max_alphabet_num - min_boundary_marker + 1) / 2;
    }
    return numer_of_variant_sites;
}
