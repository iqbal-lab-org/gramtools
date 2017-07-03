#include "sdsl/suffix_arrays.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <cassert>
#include "bwt_search.h"
#include <tuple>
#include <cstdint>

#include "variants.hpp"
#include "process_prg.hpp"
#include "map.hpp"
#include "skip.hpp"
#include "bidir_search_bwd.hpp"


/*
 * starts at end of the read, want to do fwd in the future
 *
 *
 * adds sa intervals (gives number of matches: sa interval gives you the position of each match in the prg via the suffix array)
 * sites -> variant markers within the sa interval: (everything between odd numbers)
 *  pair: site and allele
 * alleles: separated by even numbers above 5
 *
 *
 * std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> &sites
 * * std::pair -> one variant site
 * ** uint32_t -> variant site (the odd number character)
 * ** std::vector<int> -> each int is one allele, subset of alleles in variant site
 * (an index, starts at 1: 1 is first allele, 2 is second allele)
 *
 * * std::vector<std::pair<uint32_t, std::vector<int>>> -> close variants sites expect reads to cross over,
 * tracks order of crossed (by read) variant sites if variant sites close together
 *
 *
 * * std::list -> the list tracks each match of the read, elements of the list are a "match"
 *
 * sa_intervals <--one-to-one--> sites
*/
void bidir_search_bwd(const FM_Index &fm_index,
                      uint64_t left, uint64_t right,
                      uint64_t left_rev, uint64_t right_rev, // not used in bwd
                      const std::vector<uint8_t>::iterator fasta_pattern_begin,
                      const std::vector<uint8_t>::iterator fasta_pattern_end,
                      std::list<std::pair<uint64_t, uint64_t>> &sa_intervals,
                      std::list<std::pair<uint64_t, uint64_t>> &sa_intervals_rev, // not used in bwd
                      std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> &sites,
                      std::vector<int> &mask_a, const uint64_t maxx, bool &first_del,
                      const bool kmer_precalc_done, const VariantMarkers &variants,
                      std::unordered_map<uint8_t,vector<uint64_t>>& rank_all) {

    // deals with empty (first in mapping) sa interval
    if (sa_intervals.empty()) {
        sa_intervals.push_back(std::make_pair(left, right));
        sa_intervals_rev.push_back(std::make_pair(left_rev, right_rev));

        std::vector<std::pair<uint32_t, std::vector<int>>> empty_pair_vector;
        sites.push_back(empty_pair_vector);
    }

    // allele vector for "sites" variable
    std::vector<int> allele_empty;
    allele_empty.reserve(3500);

    auto fasta_pattern_it = fasta_pattern_end;
    while (fasta_pattern_it > fasta_pattern_begin) {

        if (sa_intervals.empty())
            return;

        --fasta_pattern_it;
        uint8_t fasta_char = *fasta_pattern_it;

        auto sa_interval_it = sa_intervals.begin();
        auto sa_interval_it_rev = sa_intervals_rev.begin();
        auto sites_it = sites.begin();

        if (kmer_precalc_done or fasta_pattern_it != fasta_pattern_end - 1) {

            // loop sa interval (matches of current substring)
            auto count_sa_intervals = sa_intervals.size();
            for (auto j = 0; j < count_sa_intervals; j++) {
                process_matches_overlapping_variants(fm_index, mask_a, maxx, first_del, variants,
                                                     allele_empty, sa_interval_it, sa_interval_it_rev, sites_it,
                                                     sa_intervals, sa_intervals_rev, sites);
            }
        }

        sa_interval_it = sa_intervals.begin();
        sa_interval_it_rev = sa_intervals_rev.begin();
        sites_it = sites.begin();

        first_del = match_next_charecter(fm_index, sa_intervals, sa_intervals_rev, sites, fasta_char,
                                         sa_interval_it, sa_interval_it_rev, sites_it, first_del, rank_all);
    }
}

void process_matches_overlapping_variants(const FM_Index &fm_index, const vector<int> &mask_a, const uint64_t maxx,
                                          const bool &first_del, const VariantMarkers &variants,
                                          vector<int> &allele_empty,
                                          list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it,
                                          list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it_rev,
                                          list<std::vector<std::pair<unsigned int, std::vector<int>>>>::iterator &sites_it,
                                          list<pair<uint64_t, uint64_t>> &sa_intervals,
                                          list<pair<uint64_t, uint64_t>> &sa_intervals_rev,
                                          list<vector<pair<uint32_t, vector<int>>>> &sites) {

    // check for edge of variant site
    auto sa_interval_start = (*sa_interval_it).first;
    auto sa_interval_end = (*sa_interval_it).second - 1;
    MarkerPositions marker_positions(sa_interval_start, sa_interval_end, fm_index, variants);

    uint64_t previous_marker = 0;
    uint64_t last_begin = 0;
    bool second_to_last = false;

    // loop through variant site edges (even and odd markers)
    for (auto marker_position = marker_positions.begin();
         marker_position != marker_positions.end();
         ++marker_position) {

        uint64_t marker_idx;
        uint64_t marker;
        uint64_t left_new;
        uint64_t right_new;
        uint64_t left_rev_new;
        uint64_t right_rev_new;
        bool ignore;

        add_sa_interval_for_skip(previous_marker, sa_interval_it, sa_interval_it_rev, last_begin,
                                 second_to_last, marker_position,
                                 marker_idx, marker, left_new, right_new, left_rev_new,
                                 right_rev_new, ignore);

        // takes all suffix at edge of variant and adds variant charecter to them
        // ac6cc6at5agt -> 5ac6cc6at5agt
        // last -> is end of variant site marker or not
        bool last = skip(fm_index, left_new, right_new,
                         left_rev_new, right_rev_new, marker, maxx);

        update_sites_crossed_by_reads(fm_index, sa_intervals, sa_intervals_rev, sites, mask_a,
                                      first_del, allele_empty, second_to_last, marker_idx,
                                      marker, left_new, right_new, left_rev_new, right_rev_new, ignore,
                                      last, sa_interval_it, sa_interval_it_rev, sites_it, last_begin, marker_position);

        previous_marker = marker;
    }
    ++sa_interval_it;
    ++sa_interval_it_rev;
    ++sites_it;
}


void add_sa_interval_for_skip(uint64_t previous_marker,
                              list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it,
                              list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it_rev,
                              uint64_t &last_begin, bool &second_to_last,
                              MarkerPositions &marker_it,
                              uint64_t &marker_idx, uint64_t &marker, uint64_t &left_new, uint64_t &right_new,
                              uint64_t &left_rev_new, uint64_t &right_rev_new, bool &ignore) {

    marker_idx = (*marker_it).first;
    marker = (*marker_it).second;
    left_new = (*sa_interval_it).first;
    right_new = (*sa_interval_it).second;
    left_rev_new = (*sa_interval_it_rev).first;
    right_rev_new = (*sa_interval_it_rev).second;
    ignore = ((marker == previous_marker) && (marker % 2 == 0)) ||
             ((marker % 2 == 0) && (marker == previous_marker + 1) && (marker == last_begin + 1));

    // if marker start or end of variant region
    if ((marker % 2 == 1) && (marker != previous_marker)) {
        second_to_last = false;
        last_begin = 0;
    }
    if (marker % 2 == 1) {
        left_new = marker_idx;
        right_new = marker_idx + 1;
    }
}


bool match_next_charecter(const FM_Index &fm_index, list<pair<uint64_t, uint64_t>> &sa_intervals,
                          list<pair<uint64_t, uint64_t>> &sa_intervals_rev,
                          list<vector<pair<uint32_t, vector<int>>>> &sites,
                          uint8_t fasta_char,
                          list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it,
                          list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it_rev,
                          list<std::vector<std::pair<unsigned int, std::vector<int>>>>::iterator &sites_it,
                          bool first_del, std::unordered_map<uint8_t,vector<uint64_t>>& rank_all) {

    // adds next charecter in the read; deletes sa intervals which dont match the new next character
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
    return first_del;
}


// TODO: rename to get_variant_site_edge?
std::pair<uint32_t, std::vector<int>> get_location(const FM_Index &fm_index,
                                                   const uint64_t marker_idx, const uint64_t marker,
                                                   const bool last, std::vector<int> &allele,
                                                   const std::vector<int> &mask_a) {
    uint64_t site_edge_marker;

    bool marker_is_site_edge = marker % 2 == 1;
    if (marker_is_site_edge) {
        site_edge_marker = marker;
        if (!last)
            allele.push_back(1);
    } else {
        site_edge_marker = marker - 1;
        allele.push_back(mask_a[fm_index[marker_idx]]);
    }
    return std::make_pair(site_edge_marker, allele);
}


void update_sites_crossed_by_reads(const FM_Index &fm_index, list<pair<uint64_t, uint64_t>> &sa_intervals,
                                   list<pair<uint64_t, uint64_t>> &sa_intervals_rev,
                                   list<vector<pair<uint32_t, vector<int>>>> &sites, const vector<int> &mask_a,
                                   const bool &first_del, vector<int> &allele_empty, bool &second_to_last,
                                   uint64_t marker_idx, uint64_t marker, uint64_t left_new, uint64_t right_new,
                                   uint64_t left_rev_new, uint64_t right_rev_new, bool ignore, bool last,
                                   list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it,
                                   list<std::pair<unsigned long, unsigned long>>::iterator &sa_interval_it_rev,
                                   list<std::vector<std::pair<unsigned int, std::vector<int>>>>::iterator &sites_it,
                                   uint64_t &last_begin, MarkerPositions &marker_it) {

    if (!last && (marker % 2 == 1)) {
        last_begin = marker;
        second_to_last = ((marker_it.is_second_to_last())
                          && (marker == marker_it.next_position.second));
    }

    if (sa_interval_it == sa_intervals.begin() && first_del == false && !ignore) {
        auto sa_interval = make_pair(left_new, right_new);
        sa_intervals.push_back(sa_interval);
        auto sa_interval_rev = make_pair(left_rev_new, right_rev_new);
        sa_intervals_rev.push_back(sa_interval_rev);

        auto location = get_location(fm_index, marker_idx, marker,
                                     last, allele_empty, mask_a);
        vector<pair<uint32_t, vector<int>>> tmp(1, location);
        sites.push_back(tmp);
        sites.back().reserve(100);

        allele_empty.clear();
        allele_empty.reserve(3500);

        return;
    }

    // there will be entries with pair.second empty (corresp to allele)
    // coming from crossing the last marker
    // can delete them here or in top a fcn when calculating coverages
    if (ignore) {
        if ((marker == last_begin + 1) && second_to_last) {
            auto vec_item = *(prev(sites.end(), 2));
            vec_item.back() = get_location(fm_index, marker_idx, marker, last,
                                           vec_item.back().second, mask_a);
        } else
            sites.back().back() = get_location(fm_index, marker_idx, marker, last,
                                               sites.back().back().second, mask_a);
        return;
    }

    *sa_interval_it = make_pair(left_new, right_new);
    *sa_interval_it_rev = make_pair(left_rev_new, right_rev_new);

    if ((*sites_it).back().first == marker || (*sites_it).back().first == marker - 1)
        (*sites_it).back() = get_location(fm_index, marker_idx, marker, last,
                                          (*sites_it).back().second, mask_a);

    else
        (*sites_it).push_back(get_location(fm_index, marker_idx, marker, last,
                                           allele_empty, mask_a));

	allele_empty.clear();
	allele_empty.reserve(3500);
}
