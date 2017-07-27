#include <sdsl/suffix_arrays.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "fm_index.hpp"
#include "ranks.hpp"


DNA_Rank calculate_ranks(const FM_Index &fm_index) {
    DNA_Rank rank_all;

    uint64_t bwt_size = fm_index.size();
    uint8_t symbols[] = {1, 2, 3, 4};

    for (auto symbol: symbols)
        rank_all[symbol - 1] = std::vector<uint64_t>(bwt_size, 0);

    std::vector<uint64_t> rank(4, 0);
    for (auto i = 0; i < bwt_size; i++) {
        auto curr_c = fm_index.bwt[i];
        if (curr_c > 0 && curr_c < 5) {
            rank[curr_c - 1] += 1;
            for (auto symbol: symbols)
                rank_all[symbol - 1][i] = rank[symbol - 1];
        } else {
            for (auto symbol: symbols)
                rank_all[symbol - 1][i] = rank_all[symbol - 1][i - 1];
        }
    }
    return rank_all;
}
