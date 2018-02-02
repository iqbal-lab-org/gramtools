#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

#include "sequence_read/seqread.hpp"
#include "common/parameters.hpp"
#include "common/utils.hpp"
#include "search/search.hpp"

#include "quasimap/coverage/types.hpp"
#include "quasimap/coverage/common.hpp"
#include "quasimap/quasimap.hpp"


QuasimapReadsStats quasimap_reads(const Parameters &parameters,
                                  const KmerIndex &kmer_index,
                                  const PRG_Info &prg_info) {
    std::cout << "Generating allele quasimap data structure" << std::endl;
    auto coverage = coverage::generate::empty_structure(prg_info);
    std::cout << "Done generating allele quasimap data structure" << std::endl;

    std::ofstream progress_file_handle(parameters.reads_progress_fpath);
    QuasimapReadsStats quasimap_stats = {};
    
    for (const auto &reads_fpath: parameters.reads_fpaths) {
        handle_read_file(progress_file_handle,
                         quasimap_stats,
                         coverage,
                         reads_fpath,
                         parameters,
                         kmer_index,
                         prg_info);
    }
    coverage::dump::all(coverage, parameters);
    return quasimap_stats;
}


void handle_read_file(std::ofstream &progress_file_handle,
                      QuasimapReadsStats &quasimap_stats,
                      Coverage &coverage,
                      const std::string &reads_fpath,
                      const Parameters &parameters,
                      const KmerIndex &kmer_index,
                      const PRG_Info &prg_info) {
    SeqRead reads(reads_fpath.c_str());

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
}


void quasimap_forward_reverse(QuasimapReadsStats &quasimap_reads_stats,
                              Coverage &coverage,
                              const Pattern &read,
                              const Parameters &parameters,
                              const KmerIndex &kmer_index,
                              const PRG_Info &prg_info) {
    bool read_mapped_exactly = quasimap_read(read, coverage, kmer_index, prg_info, parameters);
    if (read_mapped_exactly)
        ++quasimap_reads_stats.mapped_reads_count;

    auto reverse_read = reverse_compliment_read(read);
    read_mapped_exactly = quasimap_read(reverse_read, coverage, kmer_index, prg_info, parameters);
    if (read_mapped_exactly)
        ++quasimap_reads_stats.mapped_reads_count;
}


bool quasimap_read(const Pattern &read,
                   Coverage &coverage,
                   const KmerIndex &kmer_index,
                   const PRG_Info &prg_info,
                   const Parameters &parameters,
                   const uint32_t &random_seed) {
    auto kmer = get_kmer_from_read(parameters.kmers_size, read);
    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    auto read_mapped_exactly = not search_states.empty();
    if (not read_mapped_exactly)
        return read_mapped_exactly;
    auto read_length = read.size();
    coverage::record::search_states(coverage,
                                    search_states,
                                    read_length,
                                    prg_info,
                                    random_seed);
    return read_mapped_exactly;
}


Pattern get_kmer_from_read(const uint32_t &kmer_size, const Pattern &read) {
    Pattern kmer;
    auto kmer_start_it = read.begin() + read.size() - kmer_size;
    kmer.assign(kmer_start_it, read.end());
    return kmer;
}
