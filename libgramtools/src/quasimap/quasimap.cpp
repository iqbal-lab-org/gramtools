#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

#include "sequence_read/seqread.hpp"
#include "parameters.hpp"
#include "utils.hpp"
#include "search.hpp"

#include "quasimap/coverage/types.hpp"
#include "quasimap/coverage/allele_sum.hpp"
#include "quasimap/coverage/allele_base.hpp"
#include "quasimap/coverage/grouped_allele_counts.hpp"
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


void record_read_coverage(Coverage &coverage,
                          const SearchStates &search_states,
                          const uint64_t &read_length,
                          const PRG_Info &prg_info) {
    coverage::record::allele_sum(coverage, search_states);
    coverage::record::grouped_allele_counts(coverage, search_states);
    coverage::record::allele_base(coverage, search_states, read_length, prg_info);
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


Coverage generate_coverage_structure(const PRG_Info &prg_info) {
    Coverage coverage;
    coverage.allele_sum_coverage = coverage::generate::allele_sum_structure(prg_info);
    coverage.allele_base_coverage = coverage::generate::allele_base_structure(prg_info);
    coverage.grouped_allele_counts = coverage::generate::grouped_allele_counts(prg_info);
    return coverage;
}


Pattern get_kmer_from_read(const uint32_t& kmer_size, const Pattern &read) {
    Pattern kmer;
    auto kmer_start_it = read.begin() + read.size() - kmer_size;
    kmer.assign(kmer_start_it, read.end());
    return kmer;
}
