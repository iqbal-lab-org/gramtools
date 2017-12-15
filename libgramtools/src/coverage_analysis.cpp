#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

#include "sequence_read/seqread.hpp"
#include "parameters.hpp"
#include "utils.hpp"
#include "search.hpp"
#include "coverage_analysis.hpp"


QuasimapReadsStats quasimap_reads(const Parameters &parameters,
                             const KmerIndex &kmer_index,
                             const PRG_Info &prg_info) {
    std::cout << "Generating allele coverage data structure" << std::endl;
    auto allele_coverage = generate_allele_coverage_structure(prg_info);
    std::cout << "Done generating allele coverage data structure" << std::endl;

    SeqRead reads(parameters.reads_fpath.c_str());
    std::ofstream progress_file_handle(parameters.reads_progress_fpath);

    QuasimapReadsStats quasimap_stats = {0, 0, 0};

    for (const auto *const raw_read: reads) {
        if (quasimap_stats.all_reads_count % 1 == 0) {
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
                                 allele_coverage,
                                 read,
                                 parameters,
                                 kmer_index,
                                 prg_info);
    }
    dump_allele_coverage(allele_coverage, parameters);
    return quasimap_stats;
}


void quasimap_forward_reverse(QuasimapReadsStats &quasimap_reads_stats,
                              AlleleCoverage &allele_coverage,
                              const Pattern &read,
                              const Parameters &parameters,
                              const KmerIndex &kmer_index,
                              const PRG_Info &prg_info) {
    bool read_mapped_exactly = quasimap_read(read, allele_coverage,
                                             kmer_index, prg_info, parameters);
    if (read_mapped_exactly)
        ++quasimap_reads_stats.mapped_reads_count;
    
    auto reverse_read = reverse_compliment_read(read);
    read_mapped_exactly = quasimap_read(reverse_read, allele_coverage,
                                        kmer_index, prg_info, parameters);
    if (read_mapped_exactly)
        ++quasimap_reads_stats.mapped_reads_count;
}


bool quasimap_read(const Pattern &read,
                   AlleleCoverage &allele_coverage,
                   const KmerIndex &kmer_index,
                   const PRG_Info &prg_info,
                   const Parameters &parameters) {
    auto kmer = get_kmer_from_read(parameters.kmers_size, read);
    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    auto read_mapped_exactly = not search_states.empty();
    if (not read_mapped_exactly)
        return read_mapped_exactly;
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
            auto variant_site_coverage_index = (marker - min_boundary_marker) / 2;
            auto allele_coverage_index = allell_id - 1;

            allele_coverage[variant_site_coverage_index][allele_coverage_index] += 1;
        }
    }
}


void dump_allele_coverage(const AlleleCoverage &allele_coverage,
                          const Parameters &parameters) {
    std::ofstream file_handle(parameters.allele_coverage_fpath);
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
    auto min_boundary_marker = 5;
    auto numer_of_variant_sites = (prg_info.max_alphabet_num
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
