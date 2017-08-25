#include <sdsl/suffix_arrays.hpp>

#include "kmers.hpp"
#include "bwt_search.hpp"
#include "map.hpp"
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
 * Sites &sites
 * * std::pair -> one variant site
 * ** uint64_t -> variant site (the odd number character)
 * ** std::vector<int> -> each int is one allele, subset of alleles in variant site
 * (an index, starts at 1: 1 is first allele, 2 is second allele)
 *
 * * Site -> close variants sites expect reads to cross over,
 * tracks order of crossed (by read) variant sites if variant sites close together
 *
 *
 * * std::list -> the list tracks each match of the read, elements of the list are a "match"
 *
 * sa_intervals <--one-to-one--> sites
*/
void bidir_search_bwd(SA_Intervals &sa_intervals,
                      Sites &sites,
                      bool &delete_first_interval,
                      const std::vector<uint8_t>::const_iterator read_begin,
                      const std::vector<uint8_t>::const_iterator read_end,
                      const std::vector<int> &allele_mask,
                      const uint64_t max_alphabet_num,
                      const bool kmer_precalc_done,
                      const DNA_Rank &rank_all,
                      const FM_Index &fm_index) {

    auto read_it = read_end;
    while (read_it > read_begin) {
        if (sa_intervals.empty())
            return;

        --read_it;
        const uint8_t read_char = *read_it;
        assert((read_char >= 1) and (read_char <= 4));

        const bool read_char_is_last = read_it == read_end - 1;
        delete_first_interval = reduce_sa_intervals(sites,
                                                    sa_intervals,
                                                    delete_first_interval,
                                                    read_char,
                                                    read_char_is_last,
                                                    allele_mask,
                                                    max_alphabet_num,
                                                    kmer_precalc_done,
                                                    rank_all,
                                                    fm_index);
    }
}


bool reduce_sa_intervals(Sites &sites,
                         SA_Intervals &sa_intervals,
                         const bool delete_first_interval,
                         const uint8_t read_char,
                         const bool read_char_is_last,
                         const std::vector<int> &allele_mask,
                         const uint64_t max_alphabet_num,
                         const bool kmer_precalc_done,
                         const DNA_Rank &rank_all,
                         const FM_Index &fm_index) {

    bool new_delete_first_interval = delete_first_interval;

    auto sa_intervals_it = sa_intervals.begin();
    auto sites_it = sites.begin();

    if (kmer_precalc_done or !read_char_is_last) {
        // loop sa interval (matches of current substring)
        const auto count_sa_intervals = sa_intervals.size();
        for (auto j = 0; j < count_sa_intervals; j++) {

            auto &sa_interval = *sa_intervals_it;
            auto &site = *sites_it;
            process_reads_overlapping_variants(sa_intervals, sa_interval,
                                               sites, site,
                                               new_delete_first_interval,
                                               max_alphabet_num,
                                               allele_mask,
                                               fm_index);
            ++sa_intervals_it;
            ++sites_it;
        }
    }

    sa_intervals_it = sa_intervals.begin();
    sites_it = sites.begin();

    // examines next charecter in the read
    // deletes sa intervals which dont match the new next character
    while (sa_intervals_it != sa_intervals.end()
           and sites_it != sites.end()) {

        auto &sa_interval = *sa_intervals_it;
        auto reduced_sa_interval = bidir_search(read_char, sa_interval,
                                                rank_all, fm_index);

        const uint64_t reduced_sa_interval_size = reduced_sa_interval.second
                                                  - reduced_sa_interval.first;
        if (reduced_sa_interval_size == 0) {
            new_delete_first_interval = sa_interval == sa_intervals.front();

            sa_intervals_it = sa_intervals.erase(sa_intervals_it);
            sites_it = sites.erase(sites_it);
            continue;
        }

        // reduce SA interval to read_char interval
        sa_interval.first = reduced_sa_interval.first;
        sa_interval.second = reduced_sa_interval.second;

        ++sa_intervals_it;
        ++sites_it;
    }
    return new_delete_first_interval;
}


void process_reads_overlapping_variants(SA_Intervals &sa_intervals,
                                        SA_Interval &sa_interval,
                                        Sites &sites,
                                        Site &site,
                                        const bool delete_first_interval,
                                        const uint64_t max_alphabet_num,
                                        const std::vector<int> &allele_mask,
                                        const FM_Index &fm_index) {

    // check for edge of variant site
    const auto sa_interval_start = sa_interval.first;
    const auto sa_interval_end = sa_interval.second;

    auto marker_positions = fm_index.wavelet_tree.range_search_2d(sa_interval_start,
                                                                  sa_interval_end - 1,
                                                                  5, max_alphabet_num).second;

    uint64_t previous_marker = 0;
    uint64_t last_begin = 0;
    bool second_to_last = false;

    for (auto marker_position = marker_positions.begin(), zend = marker_positions.end();
         marker_position != zend;
         ++marker_position) {

        const auto marker_idx = (*marker_position).first;
        const auto marker = (*marker_position).second;

        uint64_t new_sa_start = sa_interval.first;
        uint64_t new_sa_end = sa_interval.second;
        if (marker % 2 == 1) {
            new_sa_start = marker_idx;
            new_sa_end = marker_idx + 1;
        }

        bool ignore = ((marker == previous_marker) and (marker % 2 == 0))
                      or ((marker % 2 == 0)
                          and (marker == previous_marker + 1)
                          and (marker == last_begin + 1));

        // if marker start or end of variant region
        if ((marker % 2 == 1) and (marker != previous_marker)) {
            second_to_last = false;
            last_begin = 0;
        }

        // takes all suffix at edge of variant and adds variant charecter to them
        // ac6cc6at5agt -> 5ac6cc6at5agt
        // last -> is end of variant site marker or not
        bool last = skip(new_sa_start, new_sa_end, max_alphabet_num, marker, fm_index);

        if (!last and (marker % 2 == 1)) {
            last_begin = marker;
            second_to_last = ((marker_position + 1 != zend)
                              and (marker == (*(marker_position + 1)).second));
        }

        const bool at_first_sa_interval = sa_interval == sa_intervals.front();
        bool res = update_sites_crossed_by_reads(sa_intervals, sites,
                                                 new_sa_start, new_sa_end,
                                                 at_first_sa_interval, ignore,
                                                 last, last_begin, second_to_last,
                                                 allele_mask, delete_first_interval,
                                                 marker, marker_idx, fm_index);

        if (res)
            continue;

        sa_interval = std::make_pair(new_sa_start, new_sa_end);

        VariantSiteMarker variant_site_marker = 0;
        auto &last_variant_site = site.back();

        if (!site.empty()) {
            variant_site_marker = last_variant_site.first;
        }

        if (variant_site_marker == marker or variant_site_marker == marker - 1) {
            auto &allele = last_variant_site.second;
            const auto variant_site = get_variant_site_edge(allele,
                                                            marker,
                                                            marker_idx,
                                                            allele_mask,
                                                            last,
                                                            fm_index);
            last_variant_site = variant_site;
            continue;
        }

        std::vector<int> allele_empty;
        allele_empty.reserve(3500);
        const auto variant_site = get_variant_site_edge(allele_empty,
                                                        marker,
                                                        marker_idx,
                                                        allele_mask,
                                                        last,
                                                        fm_index);
        site.push_back(variant_site);
        previous_marker = marker;
    }
}


bool update_sites_crossed_by_reads(SA_Intervals &sa_intervals,
                                   Sites &sites,
                                   const uint64_t left_new,
                                   const uint64_t right_new,
                                   const bool at_first_sa_interval,
                                   const bool ignore,
                                   const bool last,
                                   const uint64_t last_begin,
                                   const bool second_to_last,
                                   const std::vector<int> &allele_mask,
                                   const bool delete_first_interval,
                                   const uint64_t marker,
                                   const uint64_t marker_idx,
                                   const FM_Index &fm_index) {

    if (at_first_sa_interval and !delete_first_interval and !ignore) {
        auto sa_interval = std::make_pair(left_new, right_new);
        sa_intervals.push_back(sa_interval);

        std::vector<int> allele_empty;
        allele_empty.reserve(3500);
        auto variant_site = get_variant_site_edge(allele_empty,
                                                  marker,
                                                  marker_idx,
                                                  allele_mask,
                                                  last, fm_index);
        Site site(1, variant_site);
        sites.push_back(site);
        sites.back().reserve(100);

        return true;
    }

    // there will be entries with pair.second empty (corresp to allele)
    // coming from crossing the last marker
    // can delete them here or in top a fcn when calculating coverages
    if (ignore) {
        if ((marker == last_begin + 1) and second_to_last) {
            auto vec_item = *(std::prev(sites.end(), 2));
            const auto variant_site = get_variant_site_edge(vec_item.back().second,
                                                            marker,
                                                            marker_idx,
                                                            allele_mask,
                                                            last,
                                                            fm_index);
            vec_item.back() = variant_site;
        } else {
            auto &last_site = sites.back();
            auto &last_variant_site = last_site.back();
            auto &allele = last_variant_site.second;
            auto variant_site = get_variant_site_edge(allele,
                                                      marker,
                                                      marker_idx,
                                                      allele_mask,
                                                      last,
                                                      fm_index);
            last_variant_site = variant_site;
        }
        return true;
    }
    return false;
}


VariantSite get_variant_site_edge(std::vector<int> &allele,
                                  const uint64_t marker,
                                  const uint64_t marker_idx,
                                  const std::vector<int> &allele_mask,
                                  const bool last,
                                  const FM_Index &fm_index) {
    uint64_t site_edge_marker;

    bool marker_is_site_edge = marker % 2 == 1;
    if (marker_is_site_edge) {
        site_edge_marker = marker;
        if (!last)
            allele.push_back(1);
    } else {
        site_edge_marker = marker - 1;
        allele.push_back(allele_mask[fm_index[marker_idx]]);
    }
    return std::make_pair(site_edge_marker, allele);
}
