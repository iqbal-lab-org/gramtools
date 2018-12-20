#ifndef GRAMTOOLS_QUASIMAP_UTILS_HPP
#define GRAMTOOLS_QUASIMAP_UTILS_HPP

// TODO: Could this function be placed elsewhere to avoid this file entirely

namespace gram {
    /**
     * Retrieves the number of variant sites in the prg.
     */
    uint64_t get_number_of_variant_sites(const PRG_Info &prg_info);
}

#endif //GRAMTOOLS_QUASIMAP_UTILS_HPP
