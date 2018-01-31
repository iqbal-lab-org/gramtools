#ifndef GRAMTOOLS_PARAMETERS_HPP
#define GRAMTOOLS_PARAMETERS_HPP

enum class Commands {
    build,
    quasimap
};

struct Parameters {
    std::string gram_dirpath;
    std::string linear_prg_fpath;
    std::string encoded_prg_fpath;
    std::string fm_index_fpath;
    std::string site_mask_fpath;
    std::string allele_mask_fpath;
    std::string sdsl_memory_log_fpath;
    std::string kmer_index_fpath;
    uint32_t kmers_size;
    uint32_t max_read_size;

    // quasimap specific parameters
    std::vector<std::string> reads_fpaths;
    std::string reads_progress_fpath;

    std::string allele_sum_coverage_fpath;
    std::string allele_base_coverage_fpath;
    std::string grouped_allele_counts_fpath;
};

#endif //GRAMTOOLS_PARAMETERS_HPP
