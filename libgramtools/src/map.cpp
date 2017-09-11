#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

#include "parameters.hpp"
#include "bwt_search.hpp"
#include "map.hpp"
#include "bidir_search_bwd.hpp"


void update_coverage(AlleleCoverage &allele_coverage, const std::list<Site> &sites,
                     const std::list<SA_Interval> &sa_intervals, const std::vector<uint8_t> &encoded_read,
                     const int &count_char_in_variant_site,
                     const std::unordered_set<uint64_t> &repeats_variant_site_edge_markers,
                     const bool delete_first_interval, const PRG_Info &prg_info);


int quasimap_reads(AlleleCoverage &allele_coverage, const Parameters &params, KmerIndex &kmers, const PRG_Info &prg_info) {

    SeqRead reads(params.reads_fpath.c_str());

    // TODO: change to stdout, avoid file IO
    std::ofstream reads_fhandle(params.processed_reads_fpath);

    int count_mapped_reads = 0;
    int count_all_reads = 0;

    // variables which preserve some state across reads
    int count_char_in_variant_site = 0;
    std::unordered_set<uint64_t> repeats_variant_site_edge_markers;

    for (const auto *const read: reads) {
        ++count_all_reads;
        if (count_all_reads % 100000 == 0)
            reads_fhandle << count_all_reads << std::endl;

        std::cout << count_all_reads << std::endl;

        const auto encoded_read = int_encode_read(*read);
        const bool read_mapped = process_read(encoded_read, allele_coverage, count_char_in_variant_site,
                                              repeats_variant_site_edge_markers, kmers, params,
                                              prg_info);
        if (read_mapped)
            ++count_mapped_reads;
    }
    reads_fhandle.close();
    return count_mapped_reads;
}


bool process_read(const std::vector<uint8_t> &encoded_read, AlleleCoverage &allele_coverage, int &count_char_in_variant_site,
                  std::unordered_set<uint64_t> &repeats_variant_site_edge_markers, KmerIndex &kmers,
                  const Parameters &params, const PRG_Info &prg_info) {

    // TODO: avoid making this copy
    const auto kmer_part_start = encoded_read.begin() + encoded_read.size() - params.kmers_size;
    const auto kmer_part_end = encoded_read.end();
    const std::vector<uint8_t> read_kmer_part(kmer_part_start, kmer_part_end);

    const bool discard_read = discard_read_check(read_kmer_part, kmers);
    if (discard_read)
        return false;

    const auto kmer_size = params.kmers_size;
    bool read_mapped = quasimap_read(read_kmer_part, encoded_read, allele_coverage, count_char_in_variant_site,
                                     repeats_variant_site_edge_markers, kmers, kmer_size,
                                     prg_info);
    return read_mapped;
}


bool discard_read_check(const std::vector<uint8_t> &read_kmer_part, const KmerIndex &kmers) {
    bool kmer_not_in_precalc = (kmers.sa_intervals_map.find(read_kmer_part) == kmers.sa_intervals_map.end()
                                or kmers.sites_map.find(read_kmer_part) == kmers.sites_map.end());
    if (kmer_not_in_precalc)
        return true;

    const SA_Intervals &sa_intervals = kmers.sa_intervals_map.at(read_kmer_part);
    return sa_intervals.empty();
}


bool quasimap_read(const std::vector<uint8_t> &read_kmer_part, const std::vector<uint8_t> &encoded_read,
                   AlleleCoverage &allele_coverage, int &count_char_in_variant_site,
                   std::unordered_set<uint64_t> &repeats_variant_site_edge_markers, KmerIndex &kmers,
                   const int kmer_size, const PRG_Info &prg_info) {

    //kmers in ref means kmers that do not cross any variant site markers (non-DNA)
    //These are either in non-variable region, or are entirely within alleles
    //then the kmer does overlap a number, by definition.
    //no need to ignore first SA interval (if it was in the nonvar bit would ignore)
    bool kmer_is_nonvar = kmers.nonvar_kmers.find(read_kmer_part) != kmers.nonvar_kmers.end();
    bool delete_first_interval = !kmer_is_nonvar;

    auto &sites = kmers.sites_map[read_kmer_part];
    auto &sa_intervals = kmers.sa_intervals_map[read_kmer_part];

    const auto &read_begin = encoded_read.begin();
    const auto &read_end = encoded_read.begin() + encoded_read.size() - kmer_size;

    const bool kmer_precalc_done = true;
    bidir_search_bwd(sa_intervals, sites,
                     delete_first_interval,
                     kmer_precalc_done,
                     read_begin, read_end,
                     prg_info);

    const bool read_mapps_too_many_alleles = sa_intervals.size() > 100;
    if (read_mapps_too_many_alleles)
        return false;

    if (kmer_is_nonvar) {
        repeats_variant_site_edge_markers.clear();
        populate_repeats_variant_edges(repeats_variant_site_edge_markers,
                                       count_char_in_variant_site,
                                       sites,
                                       sa_intervals,
                                       prg_info);
    }

    update_coverage(allele_coverage, sites, sa_intervals, encoded_read, count_char_in_variant_site,
                    repeats_variant_site_edge_markers, delete_first_interval,
                    prg_info);
    return true;
}


void populate_repeats_variant_edges(std::unordered_set<uint64_t> &repeats_variant_site_edge_markers,
                                    int &count_char_in_variant_site,
                                    const Sites &sites,
                                    const SA_Intervals &sa_intervals,
                                    const PRG_Info &prg_info) {

    //becasue matches are all in non variable part of PRG
    const auto &first_site = sites.front();
    assert(first_site.empty());

    count_char_in_variant_site = 0;
    const auto &first_sa_interval = sa_intervals.front();
    for (auto ind = first_sa_interval.first; ind < first_sa_interval.second; ind++) {
        const auto charecter = prg_info.fm_index[ind];

        bool charecter_is_outside_allele = prg_info.allele_mask[charecter] == 0;
        if (charecter_is_outside_allele)
            continue;

        count_char_in_variant_site++;
        const auto variant_site_edge_marker = prg_info.sites_mask[charecter];
        repeats_variant_site_edge_markers.insert(variant_site_edge_marker);
    }
}


void update_coverage_from_sa_interval(AlleleCoverage &allele_coverage, const SA_Interval &sa_interval,
                                      const uint64_t sa_interval_size, const int &count_char_in_variant_site,
                                      const uint64_t count_repeats_variant_site_edges,
                                      const uint64_t total_num_sa_intervals, const PRG_Info &prg_info) {

    for (auto ind = sa_interval.first; ind < sa_interval.second; ind++) {
        const auto &k = prg_info.fm_index[ind];
        const auto allele = prg_info.allele_mask[k];
        if (allele == 0)
            continue;

        const auto denominator = sa_interval_size
                                 - count_char_in_variant_site
                                 + count_repeats_variant_site_edges
                                 + total_num_sa_intervals
                                 - 1;
        assert(denominator > 0);
        const auto variant_marker_coverage_idx = (allele - 5) / 2;
        auto &coverage = allele_coverage[variant_marker_coverage_idx][allele - 1];
        coverage = (coverage + 1.0) / denominator;
    }
}



void update_coverage(AlleleCoverage &allele_coverage, const std::list<Site> &sites,
                     const std::list<SA_Interval> &sa_intervals, const std::vector<uint8_t> &encoded_read,
                     const int &count_char_in_variant_site,
                     const std::unordered_set<uint64_t> &repeats_variant_site_edge_markers,
                     const bool delete_first_interval, const PRG_Info &prg_info) {

    const uint64_t count_repeats_variant_site_edges = repeats_variant_site_edge_markers.size();

    const auto &first_site = sites.front();
    const bool first_site_empty = first_site.empty();
    const uint64_t total_num_sa_intervals = sa_intervals.size();

    auto sites_it = sites.begin();
    auto sa_intervals_it = sa_intervals.begin();
    while (sa_intervals_it != sa_intervals.end() and sites_it != sites.end()) {

        const auto &site = *sites_it;
        const auto &sa_interval = *sa_intervals_it;
        const bool sa_interval_is_first = sa_interval == sa_intervals.front();
        const uint64_t sa_interval_size = sa_interval.second - sa_interval.first;

        if (sa_interval_is_first and first_site_empty) {
            //assert(!delete_first_interval); // kmer_found_in_precalc is true
            update_coverage_from_sa_interval(allele_coverage, sa_interval, sa_interval_size, count_char_in_variant_site,
                                             count_repeats_variant_site_edges, total_num_sa_intervals, prg_info);
            continue;
        }

        /*
        if (sa_interval_is_first or delete_first_interval) {
            ++sa_intervals_it;
            ++sites_it;
            continue;
        }
         */

        update_site_sa_interval_coverage(allele_coverage, site, first_site_empty, sa_interval, sa_interval_is_first,
                                         delete_first_interval, total_num_sa_intervals, count_char_in_variant_site,
                                         encoded_read, count_repeats_variant_site_edges,
                                         prg_info);

        ++sa_intervals_it;
        ++sites_it;
    }
}


void update_site_sa_interval_coverage(AlleleCoverage &allele_coverage, const Site &site, const bool first_site_empty,
                                      const SA_Interval &sa_interval, const bool sa_interval_is_first,
                                      const bool delete_first_interval, const uint64_t total_num_sa_intervals,
                                      const int count_char_in_variant_site, const std::vector<uint8_t> &encoded_read,
                                      const uint64_t count_repeats_variant_site_edges, const PRG_Info &prg_info) {

    const uint64_t sa_interval_size = sa_interval.second - sa_interval.first;

    for (const auto &variant_site: site) {
        auto variant_site_marker = variant_site.first;
        const auto &alleles = variant_site.second;

        if (variant_site != site.back() and alleles.empty())
            return;

        if (variant_site == site.back() and alleles.empty()) {
            for (auto ind = sa_interval.first; ind < sa_interval.second; ind++) {
                const auto k = prg_info.fm_index[ind];
                const auto &variant_allele_idx = prg_info.allele_mask[k] - 1;
                if (variant_allele_idx <= 0)
                    continue;

                uint64_t denominator = 0;
                if (delete_first_interval) {
                    denominator = total_num_sa_intervals;
                } else {
                    denominator = total_num_sa_intervals
                                  + sa_interval_size
                                  + count_repeats_variant_site_edges
                                  - count_char_in_variant_site
                                  - 1;
                }
                assert(denominator > 0);
                const auto variant_marker_coverage_idx = (variant_site_marker - 5) / 2;
                auto &coverage = allele_coverage[variant_marker_coverage_idx][variant_allele_idx];
                coverage = (coverage + 1.0) / denominator;
            }
            continue;
        }

        for (const auto &allele: alleles) {
            const auto variant_allele_idx = allele - 1;
            uint64_t denominator = 0;
            if (delete_first_interval) {
                denominator = total_num_sa_intervals;
            } else {
                denominator = total_num_sa_intervals
                              + sa_interval_size
                              + count_repeats_variant_site_edges
                              - count_char_in_variant_site
                              - 1;
            }
            assert(denominator > 0);
            const auto variant_marker_coverage_idx = (variant_site_marker - 5) / 2;
            auto &coverage = allele_coverage[variant_marker_coverage_idx][variant_allele_idx];
            coverage = (coverage + 1.0) / denominator;
        }
    }
}


std::vector<uint8_t> int_encode_read(const GenomicRead &read_sequence) {
    std::vector<uint8_t> encoded_read;
    encoded_read.reserve(200);

    const auto sequence_length = strlen(read_sequence.seq);

    // TODO: change to switch
    for (int i = 0; i < sequence_length; i++) {
        if (read_sequence.seq[i] == 'A' or read_sequence.seq[i] == 'a')
            encoded_read.push_back(1);
        else if (read_sequence.seq[i] == 'C' or read_sequence.seq[i] == 'c')
            encoded_read.push_back(2);
        else if (read_sequence.seq[i] == 'G' or read_sequence.seq[i] == 'g')
            encoded_read.push_back(3);
        else if (read_sequence.seq[i] == 'T' or read_sequence.seq[i] == 't')
            encoded_read.push_back(4);
        else
            break;
    }
    return encoded_read;
}


void output_allele_coverage(Parameters &params, AlleleCoverage &allele_coverage) {
    std::ofstream allele_coverage_fhandle(params.allele_coverage_fpath);
    for (uint64_t i = 0; i < allele_coverage.size(); i++) {
        for (uint64_t j = 0; j < allele_coverage[i].size(); j++)
            allele_coverage_fhandle << allele_coverage[i][j] << " ";
        allele_coverage_fhandle << std::endl;
    }
}
