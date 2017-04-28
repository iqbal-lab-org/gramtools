#ifndef GRAMTOOLS_MAP_H
#define GRAMTOOLS_MAP_H

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


Parameters parse_command_line_parameters(int argc, const char *const *argv);
void timestamp();

#endif //GRAMTOOLS_MAP_H
