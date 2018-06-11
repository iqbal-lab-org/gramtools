#include <fstream>

#include "search/search.hpp"
#include "quasimap/utils.hpp"


using namespace gram;


uint64_t gram::get_number_of_variant_sites(const PRG_Info &prg_info) {
    auto min_boundary_marker = 5;
    uint64_t number_of_variant_sites;
    if (prg_info.max_alphabet_num <= 4) {
        number_of_variant_sites = 0;
    } else {
        number_of_variant_sites = (prg_info.max_alphabet_num - min_boundary_marker + 1) / 2;
    }
    return number_of_variant_sites;
}
