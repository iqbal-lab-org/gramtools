#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>
#include <omp.h>

#include "sequence_read/seqread.hpp"

#include "common/timer_report.hpp"
#include "common/parameters.hpp"
#include "common/utils.hpp"

#include "search/search.hpp"

#include "quasimap/coverage/types.hpp"
#include "quasimap/coverage/common.hpp"
#include "quasimap/quasimap.hpp"
#include "kmer_index/load.hpp"


using namespace gram;


void commands::quasimap::run(const Parameters &parameters) {
    std::cout << "Executing quasimap command" << std::endl;
    auto timer = TimerReport();

    timer.start("Load data");
    std::cout << "Loading PRG data" << std::endl;
    const auto prg_info = load_prg_info(parameters);
    std::cout << "Loading kmer index data" << std::endl;
    const auto kmer_index = kmer_index::load(parameters);
    timer.stop();

    std::cout << "Running quasimap" << std::endl;
    timer.start("Quasimap");
    auto quasimap_stats = quasimap_reads(parameters, kmer_index, prg_info);

    std::cout << std::endl;
    std::cout << "The following counts include generated reverse complement reads."
              << std::endl;
    std::cout << "Count all reads: " << quasimap_stats.all_reads_count << std::endl;
    std::cout << "Count skipped reads: " << quasimap_stats.skipped_reads_count << std::endl;
    std::cout << "Count mapped reads: " << quasimap_stats.mapped_reads_count << std::endl;
    timer.stop();

    timer.report();
}

QuasimapReadsStats gram::quasimap_reads(const Parameters &parameters,
                                        const KmerIndex &kmer_index,
                                        const PRG_Info &prg_info) {
    std::cout << "Generating allele quasimap data structure" << std::endl;
    // The coverage structure records mapped allele counts (per site), aggregated from all mapped reads
    auto coverage = coverage::generate::empty_structure(prg_info);
    std::cout << "Done generating allele quasimap data structure" << std::endl;

    std::cout << "Processing reads:" << std::endl;
    // QuasimapReadsStats records counts of processed reads, skipped reads and mapped reads
    QuasimapReadsStats quasimap_stats = {};

    // Execute quasimap for each read file provided
    for (const auto &reads_fpath: parameters.reads_fpaths) {
        handle_read_file(quasimap_stats,
                         coverage,
                         reads_fpath,
                         parameters,
                         kmer_index,
                         prg_info);
    }
    // Write coverage results to disk
    coverage::dump::all(coverage, parameters);
    return quasimap_stats;
}


/**
 * Emplace up to `max_set_size` reads into the reads buffer.
 * Returns a vector of `Pattern`s: a `Pattern` being a vector of `Base`s, which are integer encoded.
 * The encoding of DNA letters to integers also performed in this function.
 */
std::vector<Pattern> get_reads_buffer(SeqRead::SeqIterator &reads_it, SeqRead &reads, const uint64_t &max_set_size) {
    std::vector<Pattern> reads_buffer;
    while (reads_it != reads.end() and reads_buffer.size() < max_set_size) {
        const auto *const raw_read = *reads_it;
        auto read = encode_dna_bases(*raw_read);
        reads_buffer.emplace_back(read);
        ++reads_it;
    }
    return reads_buffer;
}

/**
 * Calls the (forward_reverse) mapping routine for each read in the read buffer, in parallel (if the CL option has been specified).
 */
void handle_reads_buffer(QuasimapReadsStats &quasimap_stats,
                         Coverage &coverage,
                         const std::vector<Pattern> &reads_buffer,
                         const Parameters &parameters,
                         const KmerIndex &kmer_index,
                         const PRG_Info &prg_info) {
    uint64_t last_count_reported = 0;

    //  Parallelise loop below
    #pragma omp parallel for
    for (int i = 0; i < reads_buffer.size(); ++i) {

        auto thread_id = omp_get_thread_num();
        //  Report total number of mapped reads everytime at least `diff` such have been mapped
        if (thread_id == 0) {
            uint64_t diff = quasimap_stats.all_reads_count - last_count_reported;
            if (diff >= 10000) {
                std::cout << quasimap_stats.all_reads_count << std::endl;
                last_count_reported = quasimap_stats.all_reads_count;
            }
        }

        //  atomic: for manipulating a static variable (shared among the threads)
        #pragma omp atomic
        quasimap_stats.all_reads_count += 2; //  Increment by 2: mapping forward and reverse of read

        auto read = reads_buffer[i];
        if (read.empty()) {
            #pragma omp atomic
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

void gram::handle_read_file(QuasimapReadsStats &quasimap_stats,
                            Coverage &coverage,
                            const std::string &reads_fpath,
                            const Parameters &parameters,
                            const KmerIndex &kmer_index,
                            const PRG_Info &prg_info) {
    //  Number of reads to load in memory; is upper limit of number of reads that can be mapped in parallel
    uint64_t max_set_size = 5000;
    SeqRead reads(reads_fpath.c_str());
    auto reads_it = reads.begin();
    while (reads_it != reads.end()) {
        auto reads_buffer = get_reads_buffer(reads_it, reads, max_set_size);
        handle_reads_buffer(quasimap_stats,
                            coverage,
                            reads_buffer,
                            parameters,
                            kmer_index,
                            prg_info);
    }
}

void gram::quasimap_forward_reverse(QuasimapReadsStats &quasimap_reads_stats,
                                    Coverage &coverage,
                                    const Pattern &read,
                                    const Parameters &parameters,
                                    const KmerIndex &kmer_index,
                                    const PRG_Info &prg_info) {
    // Forward mapping
    bool read_mapped_exactly = quasimap_read(read, coverage, kmer_index, prg_info, parameters);
    if (read_mapped_exactly) {
        #pragma omp atomic
        ++quasimap_reads_stats.mapped_reads_count;
    }

    auto reverse_read = reverse_complement_read(read);
    // Reverse mapping
    read_mapped_exactly = quasimap_read(reverse_read, coverage, kmer_index, prg_info, parameters);
    if (read_mapped_exactly) {
        #pragma omp atomic
        ++quasimap_reads_stats.mapped_reads_count;
    }
}

bool gram::quasimap_read(const Pattern &read,
                         Coverage &coverage,
                         const KmerIndex &kmer_index,
                         const PRG_Info &prg_info,
                         const Parameters &parameters,
                         const uint32_t &random_seed) {
    auto kmer = get_kmer_from_read(parameters.kmers_size, read); // Gets last k bases of read

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    auto read_mapped_exactly = not search_states.empty();
    // Test read did not map
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


Pattern gram::get_kmer_from_read(const uint32_t &kmer_size, const Pattern &read) {
    Pattern kmer;
    auto kmer_start_it = read.begin() + read.size() - kmer_size;
    kmer.assign(kmer_start_it, read.end());
    return kmer;
}
