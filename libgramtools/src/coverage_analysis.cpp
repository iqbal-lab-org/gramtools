#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

#include "sequence_read/seqread.hpp"
#include "parameters.hpp"
#include "utils.hpp"
#include "search.hpp"
#include "coverage_analysis.hpp"


QuasimapStats quasimap_reads(const Parameters &params,
                             const KmerIndex &kmer_index,
                             const PRG_Info &prg_info) {
    std::cout << "Generating allele coverage data structure" << std::endl;
    auto allele_coverage = generate_allele_coverage_structure(prg_info);
    std::cout << "Done generating allele coverage data structure" << std::endl;

    SeqRead reads(params.reads_fpath.c_str());
    std::ofstream progress_file_handle(params.reads_progress_fpath);

    uint64_t all_reads_count = 0;
    uint64_t skipped_reads_count = 0;
    uint64_t mapped_reads_count = 0;

    for (const auto *const raw_read: reads) {
        if (all_reads_count % 1 == 0) {
            progress_file_handle << all_reads_count << std::endl;
            std::cout << "Reads processed: "
                      << all_reads_count
                      << std::endl;
        }
        all_reads_count++;

        auto read = encode_dna_bases(*raw_read);
        if (read.empty()) {
            ++skipped_reads_count;
            continue;
        }
        bool read_mappred_exactly = quasimap_read(read,
                                                  allele_coverage,
                                                  kmer_index,
                                                  prg_info,
                                                  params);
        if (read_mappred_exactly)
            ++mapped_reads_count;
    }
    dump_allele_coverage(allele_coverage, params);
    return std::make_tuple(all_reads_count,
                           skipped_reads_count,
                           mapped_reads_count);
}


bool quasimap_read(const Pattern &read,
                   AlleleCoverage &allele_coverage,
                   const KmerIndex &kmer_index,
                   const PRG_Info &prg_info,
                   const Parameters &params) {
    const auto kmer = get_kmer_from_read(params.kmers_size, read);
    const auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    const bool read_mapped_exactly = not search_states.empty();
    if (read_mapped_exactly)
        record_read_coverage(allele_coverage, search_states);
    return read_mapped_exactly;
}


void record_read_coverage(AlleleCoverage &allele_coverage,
                          const SearchStates &search_states) {
    for (const auto &search_state: search_states) {
        for (const auto &variant_site: search_state.variant_site_path) {
            auto marker = variant_site.first;
            auto allell_id = variant_site.second;

            auto min_boundary_marker = 5;
            auto variant_site_cover_index = (marker - min_boundary_marker) / 2;
            auto allele_cover_index = allell_id - 1;

            allele_coverage[variant_site_cover_index][allele_cover_index] += 1;
        }
    }
}


void dump_allele_coverage(const AlleleCoverage &allele_coverage,
                          const Parameters &params) {
    std::ofstream file_handle(params.allele_coverage_fpath);
    for (const auto &variant_site_coverage: allele_coverage) {
        auto allele_count = 0;
        for (const auto &coverage: variant_site_coverage) {
            file_handle << coverage;

            auto not_last_coverage = allele_count++ < variant_site_coverage.size() - 1;
            if (not_last_coverage)
                file_handle << " ";
        }
        file_handle << std::endl;
    }
}


AlleleCoverage generate_allele_coverage_structure(const PRG_Info &prg_info) {
    const auto min_boundary_marker = 5;
    const auto numer_of_variant_sites = (prg_info.max_alphabet_num
                                         - min_boundary_marker
                                         + 1)
                                        / 2;

    AlleleCoverage allele_coverage(numer_of_variant_sites);
    bool last_char_was_zero = true;

    for (const auto &prg_char: prg_info.sites_mask) {

        if (prg_char == 0) {
            last_char_was_zero = true;
            continue;
        }

        const auto &current_marker = prg_char;
        if (last_char_was_zero) {
            auto variant_site_cover_index = (current_marker - min_boundary_marker) / 2;
            allele_coverage[variant_site_cover_index].push_back(0);

            last_char_was_zero = false;
        }
    }
    return allele_coverage;
}


Pattern get_kmer_from_read(const uint32_t kmer_size, const Pattern &read) {
    Pattern kmer;
    auto kmer_start_it = read.begin() + read.size() - kmer_size;
    kmer.assign(kmer_start_it, read.end());
    return kmer;
}
