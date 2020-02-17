/** @file
 * Defines gramtools back-end commands and prg-related filepaths.
 */
#include <string>
#include <vector>


#ifndef GRAMTOOLS_PARAMETERS_HPP
#define GRAMTOOLS_PARAMETERS_HPP

namespace gram {

    enum class Commands {
        build,
        genotype
    };

    enum class Ploidy{Haploid, Diploid};

    // The number of bytes to use for each integer going on disk representing PRG string `Marker`s
    constexpr uint8_t num_bytes_per_integer{4};
    /**
     * PRG file path parameters.
     * Used for either:
     * * serialising all the necessary information for vBWT mapping to a given prg after `build` process.
     * * loading such information for `quasimap`ping reads to the prg.
     * @see gram::commands::build::parse_parameters()
     * @see gram::command::quasimap::parse_parameters()
     */
    struct Parameters {
        std::string gram_dirpath;
        std::string encoded_prg_fpath;
        std::string fm_index_fpath;
        std::string cov_graph_fpath;
        std::string sites_mask_fpath;
        std::string allele_mask_fpath;
        std::string sdsl_memory_log_fpath;

        // kmer index file paths
        std::string kmer_index_fpath;
        std::string kmers_fpath;
        std::string kmers_stats_fpath;
        std::string sa_intervals_fpath;
        std::string paths_fpath;

        Ploidy ploidy;
        uint32_t kmers_size;
        uint32_t max_read_size;
        bool all_kmers_flag;

        // quasimap specific parameters
        std::vector<std::string> reads_fpaths;

        std::string allele_sum_coverage_fpath;
        std::string allele_base_coverage_fpath;
        std::string grouped_allele_counts_fpath;

        std::string genotyped_json;
        std::string personalised_reference;

        std::string read_stats_fpath;

        uint32_t maximum_threads;
        uint32_t seed;
    };

}

#endif //GRAMTOOLS_PARAMETERS_HPP
