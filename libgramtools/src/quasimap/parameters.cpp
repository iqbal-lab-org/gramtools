#include "common/utils.hpp"
#include "quasimap/quasimap.hpp"
#include "quasimap/parameters.hpp"


Parameters commands::quasimap::parse_parameters(po::variables_map &vm,
                                                const po::parsed_options &parsed) {
    po::options_description quasimap_description("quasimap options");
    quasimap_description.add_options()
                                ("gram", po::value<std::string>(),
                                 "gramtools directory")
                                ("reads", po::value<std::vector<std::string>>()->multitoken(),
                                 "file contining reads (FASTA or FASTQ)")
                                ("kmer-size", po::value<uint32_t>(),
                                 "kmer size used in constructing the kmer index")
                                ("run-directory", po::value<std::string>(),
                                 "a directory which contains all quasimap output files")
                                ("max-threads", po::value<uint32_t>()->default_value(1),
                                 "maximum number of threads used");

    std::vector<std::string> opts = po::collect_unrecognized(parsed.options,
                                                             po::include_positional);
    opts.erase(opts.begin());
    po::store(po::command_line_parser(opts).options(quasimap_description).run(), vm);

    std::string gram_dirpath = vm["gram"].as<std::string>();

    Parameters parameters = {};
    parameters.gram_dirpath = gram_dirpath;
    parameters.linear_prg_fpath = full_path(gram_dirpath, "prg");
    parameters.encoded_prg_fpath = full_path(gram_dirpath, "encoded_prg");
    parameters.fm_index_fpath = full_path(gram_dirpath, "fm_index");
    parameters.site_mask_fpath = full_path(gram_dirpath, "variant_site_mask");
    parameters.allele_mask_fpath = full_path(gram_dirpath, "allele_mask");
    parameters.kmer_index_fpath = full_path(gram_dirpath, "kmer_index");
    parameters.kmers_fpath = full_path(gram_dirpath, "kmers");
    parameters.kmers_stats_fpath = full_path(gram_dirpath, "kmers_stats");
    parameters.sa_intervals_fpath = full_path(gram_dirpath, "sa_intervals");
    parameters.paths_fpath = full_path(gram_dirpath, "paths");

    parameters.kmers_size = vm["kmer-size"].as<uint32_t>();
    parameters.reads_fpaths = vm["reads"].as<std::vector<std::string>>();

    std::string run_dirpath = vm["run-directory"].as<std::string>();

    parameters.reads_progress_fpath = full_path(run_dirpath, "reads_progress");
    parameters.sdsl_memory_log_fpath = full_path(run_dirpath, "sdsl_memory_log");

    parameters.allele_sum_coverage_fpath = full_path(run_dirpath, "allele_sum_coverage");
    parameters.allele_base_coverage_fpath = full_path(run_dirpath, "allele_base_coverage.json");
    parameters.grouped_allele_counts_fpath = full_path(run_dirpath, "grouped_allele_counts_coverage.json");

    parameters.maximum_threads = vm["max-threads"].as<uint32_t>();
    return parameters;
}