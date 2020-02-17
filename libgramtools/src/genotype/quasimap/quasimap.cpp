#include "genotype/quasimap/quasimap.hpp"

#include "genotype/quasimap/search/BWT_search.hpp"
#include "genotype/quasimap/search/vBWT_jump.hpp"
#include "genotype/quasimap/coverage/common.hpp"
#include "genotype/quasimap/coverage/allele_base.hpp"

#include <omp.h>

using namespace gram;

QuasimapReadsStats gram::quasimap_reads(const GenotypeParams &parameters,
                                        const KmerIndex &kmer_index,
                                        const PRG_Info &prg_info,
                                        ReadStats &readstats) {
    QuasimapReadsStats quasimap_stats{};
    std::cout << "Generating allele quasimap data structure" << std::endl;
    // The coverage structure records mapped allele counts (per site), aggregated from all mapped reads
    quasimap_stats.coverage = coverage::generate::empty_structure(prg_info);
    std::cout << "Done generating allele quasimap data structure" << std::endl;

    std::cout << "Processing reads:" << std::endl;
    // QuasimapReadsStats records counts of processed reads, skipped reads and mapped reads

    // Execute quasimap for each read file provided
    for (const auto &reads_fpath: parameters.reads_fpaths) {
        handle_read_file(quasimap_stats,
                         reads_fpath,
                         parameters,
                         kmer_index,
                         prg_info);
    }

    auto& coverage = quasimap_stats.coverage;
    //Compute read mapping statistics (used in `infer` command). Can only be done after mapping!
    readstats.compute_coverage_depth(coverage, prg_info.coverage_graph.par_map);

    // Extract non-nested per base coverage
    coverage.allele_base_coverage = coverage::generate::allele_base_non_nested(prg_info);

    // Write coverage results to disk
    coverage::dump::all(coverage, parameters);
    return quasimap_stats;
}


/**
 * Emplace up to `max_set_size` reads into the reads buffer.
 * Returns a vector of `Pattern`s: a `Pattern` being a vector of `Base`s, which are integer encoded.
 * The encoding of DNA letters to integers also performed in this function.
 */
std::vector<Sequence> get_reads_buffer(SeqRead::SeqIterator &reads_it, SeqRead &reads, const uint64_t &max_set_size) {
    std::vector<Sequence> reads_buffer;
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
void handle_reads_buffer(QuasimapReadsStats &quasimap_stats, const std::vector<Sequence> &reads_buffer,
                         const GenotypeParams &parameters, const KmerIndex &kmer_index, const PRG_Info &prg_info) {
    uint64_t last_count_reported = 0;

    //  Parallelise loop below
    #pragma omp parallel for
    for (std::size_t i = 0; i < reads_buffer.size(); ++i) {

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
                                 read,
                                 parameters,
                                 kmer_index,
                                 prg_info);
    }
}

void
gram::handle_read_file(QuasimapReadsStats &quasimap_stats, const std::string &reads_fpath, const GenotypeParams &parameters,
                       const KmerIndex &kmer_index, const PRG_Info &prg_info) {
    //  Number of reads to load in memory; is upper limit of number of reads that can be mapped in parallel
    uint64_t max_set_size = 5000;
    SeqRead reads(reads_fpath.c_str());
    auto reads_it = reads.begin();
    while (reads_it != reads.end()) {
        auto reads_buffer = get_reads_buffer(reads_it, reads, max_set_size);
        handle_reads_buffer(quasimap_stats,
                            reads_buffer,
                            parameters,
                            kmer_index,
                            prg_info);
    }
}

void
gram::quasimap_forward_reverse(QuasimapReadsStats &quasimap_stats, const Sequence &read, const GenotypeParams &parameters,
                               const KmerIndex &kmer_index, const PRG_Info &prg_info) {
    // Forward mapping
    bool read_mapped_exactly = quasimap_read(read, quasimap_stats.coverage,
            kmer_index, prg_info, parameters);
    if (read_mapped_exactly) {
        #pragma omp atomic
        ++quasimap_stats.mapped_reads_count;
    }

    auto reverse_read = reverse_complement_read(read);
    // Reverse mapping
    read_mapped_exactly = quasimap_read(reverse_read, quasimap_stats.coverage,
            kmer_index, prg_info, parameters);
    if (read_mapped_exactly) {
        #pragma omp atomic
        ++quasimap_stats.mapped_reads_count;
    }
}

bool gram::quasimap_read(const Sequence &read,
                         Coverage &coverage,
                         const KmerIndex &kmer_index,
                         const PRG_Info &prg_info,
                         const GenotypeParams &parameters) {
    auto kmer = get_kmer_from_read(parameters.kmers_size, read); // Gets last k bases of read

    auto search_states = search_read_backwards(read, kmer, kmer_index, prg_info);
    auto read_mapped_exactly = not search_states.empty();
    // Test read did not map
    if (not read_mapped_exactly)
        return read_mapped_exactly;
    auto read_length = read.size();

    uint64_t random_seed = parameters.seed;
    coverage::record::search_states(coverage,
                                    search_states,
                                    read_length,
                                    prg_info,
                                    random_seed);
    return read_mapped_exactly;
}


Sequence gram::get_kmer_from_read(const uint32_t &kmer_size, const Sequence &read) {
    Sequence kmer;
    auto kmer_start_it = read.begin() + read.size() - kmer_size;
    kmer.assign(kmer_start_it, read.end());
    return kmer;
}

SearchStates gram::search_read_backwards(const Sequence &read,
                                         const Sequence &kmer,
                                         const KmerIndex &kmer_index,
                                         const PRG_Info &prg_info) {
    // Test if kmer has been indexed
    bool kmer_in_index = kmer_index.find(kmer) != kmer_index.end();
    if (not kmer_in_index)
        return SearchStates{};

    // Test if kmer has been indexed, but has no search states in prg
    auto kmer_index_search_states = kmer_index.at(kmer);
    if (kmer_index_search_states.empty())
        return kmer_index_search_states;


    // Reverse iterator + skipping through indexed kmer in read
    auto read_begin = read.rbegin();
    std::advance(read_begin, kmer.size());

    SearchStates new_search_states = kmer_index_search_states;

    for (auto it = read_begin; it != read.rend(); ++it) { /// Iterates end to start of read
        const int_Base &pattern_char = *it;
        new_search_states = process_read_char_search_states(pattern_char,
                                                            new_search_states,
                                                            prg_info);
        // Test if no mapping found upon character extension
        auto read_not_mapped = new_search_states.empty();
        if (read_not_mapped)
            break;
    }

    new_search_states = handle_allele_encapsulated_states(new_search_states, prg_info);
    return new_search_states;
}

SearchStates gram::process_read_char_search_states(const int_Base &pattern_char,
                                                   const SearchStates &old_search_states,
                                                   const PRG_Info &prg_info) {
    //  Before extending backward search with next character, check for variant markers in the current SA intervals
    //  This is the v part of vBWT.
    auto post_markers_search_states = process_markers_search_states(old_search_states,
                                                                    prg_info);
    //  Regular backward searching
    auto new_search_states = search_base_backwards(pattern_char,
                                                   post_markers_search_states,
                                                   prg_info);
    return new_search_states;
}


/**
 * Produce integer-encoded Watson-Crick base complement.
 */
int_Base complement_encoded_base(const int_Base &encoded_base) {
    switch (encoded_base) {
        case 1:
            return 4;
        case 2:
            return 3;
        case 3:
            return 2;
        case 4:
            return 1;
        default:
            return 0;
    }
}

Sequence gram::reverse_complement_read(const Sequence &read) {
    Sequence reverse_read;
    reverse_read.reserve(read.size());

    for (auto it = read.rbegin(); it != read.rend(); ++it) {
        const auto &base = *it;
        auto complement_base = complement_encoded_base(base);
        reverse_read.push_back(complement_base);
    }
    return reverse_read;
}
