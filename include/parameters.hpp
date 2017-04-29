#ifndef GRAMTOOLS_PARAMETERS_HPP
#define GRAMTOOLS_PARAMETERS_HPP

struct Parameters{
    std::string prg_fpath;
    std::string csa_fpath;
    std::string festa_fpath;
    std::string site_mask_fpath;
    std::string allele_mask_fpath;
    std::string allele_coverage_fpath;
    std::string processed_reads_fpath;
    std::string prg_integer_alphabet_fpath;
    std::string csa_memory_log_fpath;
    std::string prg_kmers_fpath;
    int kmers_size;
};

#endif //GRAMTOOLS_PARAMETERS_HPP
