#ifndef GRAMTOOLS_PARAMETERS_HPP
#define GRAMTOOLS_PARAMETERS_HPP

enum class Commands {
    build,
    quasimap
};

struct Parameters {
    std::string linear_prg_fpath;
    std::string encoded_prg_fpath;
    std::string fm_index_fpath;
    std::string site_mask_fpath;
    std::string allele_mask_fpath;
    std::string sdsl_memory_log_fpath;
    std::string kmer_suffix_diffs_fpath;
    std::string kmer_index_fpath;
    uint32_t kmers_size;

    // quasimap specific parameters
    std::string reads_fpath;
    std::string allele_coverage_fpath;
    std::string reads_progress_fpath;
};

#endif //GRAMTOOLS_PARAMETERS_HPP
