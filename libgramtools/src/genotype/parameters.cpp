#include "common/utils.hpp"
#include "genotype/quasimap/quasimap.hpp"
#include "genotype/parameters.hpp"


using namespace gram;


Parameters commands::genotype::parse_parameters(po::variables_map &vm,
                                                const po::parsed_options &parsed) {
    std::string gram_dirpath;
    std::string run_directory;
    uint32_t kmer_size;
    std::vector<std::string> reads;

    po::options_description quasimap_description("genotype options");
    quasimap_description.add_options()
                                ("gram", po::value<std::string>(&gram_dirpath)->required(),
                                 "gramtools directory")
                                ("reads", po::value<std::vector<std::string>>(&reads)->multitoken()->required(),
                                 "file containing reads (FASTA or FASTQ)")
                                ("kmer-size", po::value<uint32_t>(&kmer_size)->required(),
                                 "kmer size used in constructing the kmer index")
                                ("run-directory", po::value<std::string>(&run_directory)->required(),
                                 "the directory where to store all quasimap output files")
                                ("max-threads", po::value<uint32_t>()->default_value(1),
                                 "maximum number of threads used")
                                ("seed", po::value<uint32_t>()->default_value(0),
                                        "seed for pseudo-random selection of multi-mapping reads. the default of 0 produces a random seed.");

    std::vector<std::string> opts = po::collect_unrecognized(parsed.options,
                                                             po::include_positional);
    if (opts.size() > 0) opts.erase(opts.begin()); // Takes out the command itself
    try {
        po::store(po::command_line_parser(opts).options(quasimap_description).run(), vm);
        po::notify(vm);
    }
    catch(const std::exception &e){
        std::cout << e.what() << std::endl;
        std::cout << quasimap_description << std::endl;
        exit(1);
    }

    Parameters parameters = {};
    parameters.gram_dirpath = gram_dirpath;
    parameters.encoded_prg_fpath = full_path(gram_dirpath, "prg");
    parameters.fm_index_fpath = full_path(gram_dirpath, "fm_index");
    parameters.cov_graph_fpath = full_path(gram_dirpath, "cov_graph");
    parameters.sites_mask_fpath = full_path(gram_dirpath, "variant_site_mask");
    parameters.allele_mask_fpath = full_path(gram_dirpath, "allele_mask");

    parameters.kmer_index_fpath = full_path(gram_dirpath, "kmer_index");
    parameters.kmers_fpath = full_path(gram_dirpath, "kmers");
    parameters.kmers_stats_fpath = full_path(gram_dirpath, "kmers_stats");
    parameters.sa_intervals_fpath = full_path(gram_dirpath, "sa_intervals");
    parameters.paths_fpath = full_path(gram_dirpath, "paths");

    parameters.kmers_size = kmer_size;
    parameters.reads_fpaths = reads;

    std::string run_dirpath = run_directory;
    parameters.sdsl_memory_log_fpath = full_path(run_dirpath, "sdsl_memory_log");

    parameters.allele_sum_coverage_fpath = full_path(run_dirpath, "allele_sum_coverage");
    parameters.allele_base_coverage_fpath = full_path(run_dirpath, "allele_base_coverage.json");
    parameters.grouped_allele_counts_fpath = full_path(run_dirpath, "grouped_allele_counts_coverage.json");
    
    parameters.read_stats_fpath = full_path(run_dirpath, "read_stats.json");

    parameters.maximum_threads = vm["max-threads"].as<uint32_t>();
    parameters.seed = vm["seed"].as<uint32_t>();
    return parameters;
}