#include <sdsl/suffix_arrays.hpp>

#include "bwt_search.h"
#include "map.hpp"

using namespace sdsl;


std::vector<uint8_t>::iterator bidir_search_bwd(const FM_Index &fm_index,
                                                uint64_t left, uint64_t right,
                                                uint64_t left_rev, uint64_t right_rev,
                                                const std::vector<uint8_t>::iterator fasta_pattern_begin,
                                                const std::vector<uint8_t>::iterator fasta_pattern_end,
                                                std::list<std::pair<uint64_t, uint64_t>> &sa_intervals,
                                                std::list<std::pair<uint64_t, uint64_t>> &sa_intervals_rev,
                                                std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> &sites,
                                                std::vector<int> &mask_a, const uint64_t maxx, bool &first_del,
                                                const bool kmer_precalc_done, const VariantMarkers &variants,
                                                std::unordered_map<uint8_t, vector<uint64_t>> &rank_all) {

    if (sa_intervals.empty()) {
        sa_intervals.push_back(std::make_pair(left, right));
        sa_intervals_rev.push_back(std::make_pair(left_rev, right_rev));

        std::vector<std::pair<uint32_t, std::vector<int>>> empty_pair_vector;
        sites.push_back(empty_pair_vector);
    }

    std::vector<int> allele_empty;
    allele_empty.reserve(3500);

    auto fasta_pattern_it = fasta_pattern_end;
    while (fasta_pattern_it > fasta_pattern_begin) {

        if (sa_intervals.empty())
            break;

        --fasta_pattern_it;
        uint8_t fasta_char = *fasta_pattern_it;

        auto sa_interval_it = sa_intervals.begin();
        auto sa_interval_it_rev = sa_intervals_rev.begin();
        auto sites_it = sites.begin();

        if ((fasta_pattern_it != fasta_pattern_end - 1) or kmer_precalc_done) {

            const uint64_t count_sa_intervals = sa_intervals.size();
            for (auto j = 0; j < count_sa_intervals; j++) {

                auto sa_interval_start = (*sa_interval_it).first;
                auto sa_interval_end = (*sa_interval_it).second - 1;
                auto variant_markers = find_variant_markers(sa_interval_start, sa_interval_end,
                                                            fm_index, variants);
                uint64_t previous_marker = 0;
                uint64_t last_begin = 0;
                bool second_to_last = false;

                for (auto marker_it = variant_markers.begin();
                     marker_it != variant_markers.end();
                     ++marker_it) {

                    uint64_t marker_idx = (*marker_it).first;
                    uint64_t marker = (*marker_it).second;

                    if ((marker % 2 == 1) && (marker != previous_marker)) {
                        second_to_last = false;
                        last_begin = 0;
                    }

                    auto left_new = (*sa_interval_it).first;
                    auto right_new = (*sa_interval_it).second;
                    if (marker % 2 == 1) {
                        left_new = marker_idx;
                        right_new = marker_idx + 1;
                    }

                    auto left_rev_new = (*sa_interval_it_rev).first;
                    auto right_rev_new = (*sa_interval_it_rev).second;

                    bool ignore = ((marker == previous_marker) && (marker % 2 == 0)) ||
                                  ((marker % 2 == 0) && (marker == previous_marker + 1) && (marker == last_begin + 1));

                    bool last = skip(fm_index, left_new, right_new,
                                     left_rev_new, right_rev_new, marker, maxx);
                    if (!last && (marker % 2 == 1)) {
                        last_begin = marker;
                        if ((marker_it + 1 != variant_markers.end()) && (marker == (*(marker_it + 1)).second))
                            second_to_last = true;
                    }

                    if (sa_interval_it == sa_intervals.begin() && first_del == false && !ignore) {
                        auto sa_interval = std::make_pair(left_new, right_new);
                        sa_intervals.push_back(sa_interval);
                        auto sa_interval_rev = std::make_pair(left_rev_new, right_rev_new);
                        sa_intervals_rev.push_back(sa_interval_rev);

                        auto location = get_location(fm_index, marker_idx, marker,
                                                     last, allele_empty, mask_a);
                        std::vector<std::pair<uint32_t, std::vector<int>>> tmp(1, location);
                        sites.push_back(tmp);
                        sites.back().reserve(100);

                        allele_empty.clear();
                        allele_empty.reserve(3500);
                    } else {
                        // there will be entries with pair.second empty (corresp to allele)
                        // coming from crossing the last marker
                        // can delete them here or in top a fcn when calculating coverages
                        if (ignore) {
                            if ((marker == last_begin + 1) && second_to_last) {
                                auto vec_item = *(std::prev(sites.end(), 2));
                                vec_item.back() = get_location(fm_index, marker_idx, marker, last,
                                                               vec_item.back().second, mask_a);
                            } else
                                sites.back().back() = get_location(fm_index, marker_idx, marker, last,
                                                                   sites.back().back().second, mask_a);
                        } else {
                            *sa_interval_it = std::make_pair(left_new, right_new);
                            *sa_interval_it_rev = std::make_pair(left_rev_new, right_rev_new);

                            if ((*sites_it).back().first == marker || (*sites_it).back().first == marker - 1)
                                (*sites_it).back() = get_location(fm_index, marker_idx, marker, last,
                                                                  (*sites_it).back().second, mask_a);

                            else
                                (*sites_it).push_back(get_location(fm_index, marker_idx, marker, last,
                                                                   allele_empty, mask_a));

                            allele_empty.clear();
                            allele_empty.reserve(3500);
                        }
                    }
                    previous_marker = marker;
                } // for (auto z = variant_markers.begin(), zend = variant_markers.end(); z != zend; ++z)
                ++sa_interval_it;
                ++sa_interval_it_rev;
                ++sites_it;
            } // for (uint64_t j = 0; j < count_sa_intervals; j++)
        } // if ((fasta_pattern_it != fasta_pattern_end - 1) or kmer_precalc_done)

        sa_interval_it = sa_intervals.begin();
        sa_interval_it_rev = sa_intervals_rev.begin();
        sites_it = sites.begin();

        while (sa_interval_it != sa_intervals.end()
               && sa_interval_it_rev != sa_intervals_rev.end()
               && sites_it != sites.end()) {
            //calculate sum to return- can do this in top fcns
            auto tmp2 = bidir_search(fm_index, (*sa_interval_it).first, (*sa_interval_it).second,
                                     (*sa_interval_it_rev).first, (*sa_interval_it_rev).second,
                                     fasta_char, rank_all);
            if (tmp2 > 0) {
                ++sa_interval_it;
                ++sa_interval_it_rev;
                ++sites_it;
                continue;
            }

            if (sa_interval_it == sa_intervals.begin())
                first_del = true;

            sa_interval_it = sa_intervals.erase(sa_interval_it);
            sa_interval_it_rev = sa_intervals_rev.erase(sa_interval_it_rev);
            sites_it = sites.erase(sites_it);
        }
    } // while (fasta_pattern_it > fasta_pattern_begin && !sa_intervals.empty())

    if (fasta_pattern_it != fasta_pattern_begin)
        return fasta_pattern_it;

    if (!sa_intervals.empty())
        return fasta_pattern_end;

    return fasta_pattern_begin;
}
