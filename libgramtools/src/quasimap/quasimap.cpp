#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

#include "sequence_read/seqread.hpp"
#include "parameters.hpp"
#include "utils.hpp"
#include "search.hpp"
#include "quasimap/quasimap.hpp"


QuasimapReadsStats quasimap_reads(const Parameters &parameters,
                                  const KmerIndex &kmer_index,
                                  const PRG_Info &prg_info) {
    std::cout << "Generating allele quasimap data structure" << std::endl;
    auto coverage = generate_coverage_structure(prg_info);
    std::cout << "Done generating allele quasimap data structure" << std::endl;

    SeqRead reads(parameters.reads_fpath.c_str());
    std::ofstream progress_file_handle(parameters.reads_progress_fpath);

    QuasimapReadsStats quasimap_stats;

    for (const auto *const raw_read: reads) {
        if (quasimap_stats.all_reads_count % 100 == 0) {
            progress_file_handle
                    << quasimap_stats.all_reads_count
                    << std::endl;
            std::cout << "Reads processed: "
                      << quasimap_stats.all_reads_count
                      << std::endl;
        }
        quasimap_stats.all_reads_count += 2;

        auto read = encode_dna_bases(*raw_read);
        if (read.empty()) {
            quasimap_stats.skipped_reads_count += 2;
            continue;
        }
        quasimap_forward_reverse(quasimap_stats,
                                 coverage,
                                 read,
                                 parameters,
                                 kmer_index,
                                 prg_info);
    }
    dump_coverage(coverage, parameters);
    return quasimap_stats;
}


void quasimap_forward_reverse(QuasimapReadsStats &quasimap_reads_stats,
                              Coverage &coverage,
                              const Pattern &read,
                              const Parameters &parameters,
                              const KmerIndex &kmer_index,
                              const PRG_Info &prg_info) {
    bool read_mapped_exactly = quasimap_read(read, coverage,
                                             kmer_index, prg_info, parameters);
    if (read_mapped_exactly)
        ++quasimap_reads_stats.mapped_reads_count;

    auto reverse_read = reverse_compliment_read(read);
    read_mapped_exactly = quasimap_read(reverse_read, coverage,
                                        kmer_index, prg_info, parameters);
    if (read_mapped_exactly)
        ++quasimap_reads_stats.mapped_reads_count;
}


bool quasimap_read(const Pattern &read,
                   Coverage &coverage,
                   const KmerIndex &kmer_index,
                   const PRG_Info &prg_info,
                   const Parameters &parameters) {
    auto kmer = get_kmer_from_read(parameters.kmers_size, read);
    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    auto read_mapped_exactly = not search_states.empty();
    if (not read_mapped_exactly)
        return read_mapped_exactly;
    auto read_length = read.size();
    record_read_coverage(coverage, search_states, read_length, prg_info);
    return read_mapped_exactly;
}


void record_allele_sum_coverage(Coverage &coverage,
                                const SearchStates &search_states) {
    auto &allele_sum_coverage = coverage.allele_sum_coverage;
    for (const auto &search_state: search_states) {
        for (const auto &variant_site: search_state.variant_site_path) {
            auto marker = variant_site.first;
            auto allell_id = variant_site.second;

            auto min_boundary_marker = 5;
            auto variant_site_coverage_index = (marker - min_boundary_marker) / 2;
            auto allele_coverage_index = allell_id - 1;

            allele_sum_coverage[variant_site_coverage_index][allele_coverage_index] += 1;
        }
    }
}


void record_grouped_allele_counts(Coverage &coverage,
                                  const SearchStates &search_states) {
    auto &allele_sum_coverage = coverage.allele_sum_coverage;
    for (const auto &search_state: search_states) {

        /*
        for (const auto &variant_site: search_state.variant_site_path) {
            auto marker = variant_site.first;
            auto allell_id = variant_site.second;

            auto min_boundary_marker = 5;
            auto variant_site_coverage_index = (marker - min_boundary_marker) / 2;
            auto allele_coverage_index = allell_id - 1;

            allele_sum_coverage[variant_site_coverage_index][allele_coverage_index] += 1;
        }
        */

    }
}


uint64_t set_site_base_coverage(Coverage &coverage,
                                const VariantSite &path_element,
                                const uint64_t allele_start_index,
                                const uint64_t max_set_bases) {
    auto marker = path_element.first;
    auto min_boundary_marker = 5;
    auto variant_site_coverage_index = (marker - min_boundary_marker) / 2;
    auto &site_coverage = coverage.allele_base_coverage.at(variant_site_coverage_index);

    auto allell_id = path_element.second;
    auto allele_coverage_index = allell_id - 1;
    auto &allele_coverage = site_coverage.at(allele_coverage_index);

    auto end_index = std::min(max_set_bases, allele_coverage.size());
    uint64_t count_bases_covered = end_index - allele_start_index;

    for (uint64_t i = allele_start_index; i < end_index; ++i)
        ++allele_coverage[i];
    return count_bases_covered;
}


std::pair<uint64_t, uint64_t> site_marker_prg_indexes(const uint64_t &site_marker, const PRG_Info &prg_info) {
    auto alphabet_rank = prg_info.fm_index.char2comp[site_marker];
    auto first_sa_index = prg_info.fm_index.C[alphabet_rank];
    auto second_sa_index = first_sa_index + 1;

    auto first_prg_index = prg_info.fm_index[first_sa_index];
    auto second_prg_index = prg_info.fm_index[second_sa_index];

    if (first_prg_index < second_prg_index)
        return std::make_pair(first_prg_index, second_prg_index);
    else
        return std::make_pair(second_prg_index, first_prg_index);
}


uint64_t inter_site_base_count(const uint64_t &first_site_marker,
                               const uint64_t &second_site_marker,
                               const PRG_Info &prg_info) {
    auto first_site_prg_start_end = site_marker_prg_indexes(first_site_marker, prg_info);
    auto first_site_prg_end = first_site_prg_start_end.second;

    auto second_site_prg_start_end = site_marker_prg_indexes(second_site_marker, prg_info);
    auto second_site_prg_start = second_site_prg_start_end.first;

    return second_site_prg_start - first_site_prg_end - 1;
}


uint64_t allele_start_offset_index(const uint64_t within_allele_prg_index, const PRG_Info &prg_info) {
    // find the nearest left marker with rank and select
    auto number_markers_before = prg_info.prg_markers_rank(within_allele_prg_index);
    auto marker_index = prg_info.prg_markers_select(number_markers_before);
    auto offset = within_allele_prg_index - marker_index - 1;
    return offset;
}


void sa_index_allele_base_coverage(Coverage &coverage,
                                   const uint64_t &sa_index,
                                   const uint64_t &read_length,
                                   const SearchState &search_state,
                                   const PRG_Info &prg_info) {
    uint64_t read_bases_consumed = 0;
    uint64_t last_site_marker = 0;
    auto path_it = search_state.variant_site_path.begin();

    auto read_start_index = prg_info.fm_index[sa_index];
    auto start_site_marker = prg_info.sites_mask[read_start_index];
    bool read_starts_within_site = start_site_marker != 0;
    if (read_starts_within_site) {
        const auto &path_element = *path_it;
        last_site_marker = path_element.first;

        auto allele_coverage_offset = allele_start_offset_index(read_start_index, prg_info);
        auto max_set_bases = read_length - read_bases_consumed;
        read_bases_consumed += set_site_base_coverage(coverage,
                                                      path_element,
                                                      allele_coverage_offset,
                                                      max_set_bases);

        ++path_it;
    } else {
        // forward to first path element
        const auto &path_element = *path_it;
        auto site_marker = path_element.first;
        auto site_prg_start_end = site_marker_prg_indexes(site_marker, prg_info);
        read_bases_consumed += site_prg_start_end.first - read_start_index;
    }

    auto last_path_it = search_state.variant_site_path.end();
    while (read_bases_consumed < read_length and path_it != last_path_it) {
        // take account of bases between sites
        const auto &path_element = *path_it;
        auto site_marker = path_element.first;

        if (last_site_marker != 0)
            read_bases_consumed += inter_site_base_count(last_site_marker, site_marker, prg_info);
        last_site_marker = site_marker;

        uint64_t allele_coverage_offset = 0;
        auto max_set_bases = read_length - read_bases_consumed;
        read_bases_consumed += set_site_base_coverage(coverage,
                                                      path_element,
                                                      allele_coverage_offset,
                                                      max_set_bases);
        ++path_it;
    }
}


void record_allele_base_coverage(Coverage &coverage,
                                 const SearchStates &search_states,
                                 const uint64_t &read_length,
                                 const PRG_Info &prg_info) {
    for (const auto &search_state: search_states) {
        if (search_state.variant_site_path.empty())
            continue;

        auto first_sa_index = search_state.sa_interval.first;
        auto last_sa_index = search_state.sa_interval.second;
        for (auto sa_index = first_sa_index; sa_index <= last_sa_index; ++sa_index)
            sa_index_allele_base_coverage(coverage,
                                          sa_index,
                                          read_length,
                                          search_state,
                                          prg_info);
    }
}


void record_read_coverage(Coverage &coverage,
                          const SearchStates &search_states,
                          const uint64_t &read_length,
                          const PRG_Info &prg_info) {
    record_allele_sum_coverage(coverage, search_states);
    record_grouped_allele_counts(coverage, search_states);
    record_allele_base_coverage(coverage, search_states, read_length, prg_info);
}


void dump_coverage(const Coverage &coverage,
                   const Parameters &parameters) {
    std::ofstream file_handle(parameters.allele_coverage_fpath);
    for (const auto &variant_site_coverage: coverage.allele_sum_coverage) {
        auto allele_count = 0;
        for (const auto &sum_coverage: variant_site_coverage) {
            file_handle << sum_coverage;
            auto not_last_coverage = allele_count++ < variant_site_coverage.size() - 1;
            if (not_last_coverage)
                file_handle << " ";
        }
        file_handle << std::endl;
    }
}


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


AlleleSumCoverage generate_allele_sum_coverage_structure(const PRG_Info &prg_info) {
    uint64_t numer_of_variant_sites = get_number_of_variant_sites(prg_info);
    AlleleSumCoverage allele_sum_coverage(numer_of_variant_sites);

    const auto min_boundary_marker = 5;
    bool last_char_was_zero = true;

    for (const auto &mask_value: prg_info.sites_mask) {
        if (mask_value == 0) {
            last_char_was_zero = true;
            continue;
        }

        const auto &current_marker = mask_value;
        if (last_char_was_zero) {
            auto variant_site_cover_index = (current_marker - min_boundary_marker) / 2;
            allele_sum_coverage[variant_site_cover_index].push_back(0);
            last_char_was_zero = false;
        }
    }
    return allele_sum_coverage;
}


SitesAlleleBaseCoverage generate_base_coverage_structure(const PRG_Info &prg_info) {
    uint64_t numer_of_variant_sites = get_number_of_variant_sites(prg_info);
    SitesAlleleBaseCoverage allele_base_coverage(numer_of_variant_sites);

    const auto min_boundary_marker = 5;

    uint64_t allele_size = 0;
    Marker last_marker = 0;

    for (const auto &mask_value: prg_info.sites_mask) {
        auto within_allele = mask_value != 0;
        if (within_allele) {
            allele_size += 1;
            last_marker = mask_value;
            continue;
        }

        auto no_allele_to_flush = allele_size == 0;
        if (no_allele_to_flush)
            continue;

        BaseCoverage bases(allele_size);
        uint64_t variant_site_cover_index = (last_marker - min_boundary_marker) / 2;
        allele_base_coverage.at(variant_site_cover_index).emplace_back(bases);
        allele_size = 0;
    }
    return allele_base_coverage;
}


Coverage generate_coverage_structure(const PRG_Info &prg_info) {
    Coverage coverage;
    coverage.allele_sum_coverage = generate_allele_sum_coverage_structure(prg_info);
    coverage.allele_base_coverage = generate_base_coverage_structure(prg_info);

    uint64_t numer_of_variant_sites = get_number_of_variant_sites(prg_info);
    coverage.grouped_allele_counts.reserve(numer_of_variant_sites);
    return coverage;
}


Pattern get_kmer_from_read(const uint32_t kmer_size, const Pattern &read) {
    Pattern kmer;
    auto kmer_start_it = read.begin() + read.size() - kmer_size;
    kmer.assign(kmer_start_it, read.end());
    return kmer;
}
