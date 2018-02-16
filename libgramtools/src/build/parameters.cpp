#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/variant/variant.hpp>
#include <boost/filesystem.hpp>

#include "common/parameters.hpp"
#include "common/utils.hpp"

#include "build/parameters.hpp"


Parameters commands::build::parse_parameters(po::variables_map &vm, const po::parsed_options &parsed) {
    po::options_description build_description("build options");
    build_description.add_options()
                             ("gram", po::value<std::string>(),
                              "gramtools directory")
                             ("kmer-size", po::value<uint32_t>(),
                              "kmer size used in constructing the kmer index")
                             ("max-read-size", po::value<uint32_t>(),
                              "read maximum size for the set of reads used when quasimaping");

    std::vector<std::string> opts = po::collect_unrecognized(parsed.options,
                                                             po::include_positional);
    opts.erase(opts.begin());

    po::store(po::command_line_parser(opts).options(build_description).run(), vm);

    std::string gram_dirpath = vm["gram"].as<std::string>();

    Parameters parameters;
    parameters.gram_dirpath = gram_dirpath;
    parameters.linear_prg_fpath = full_path(gram_dirpath, "prg");
    parameters.encoded_prg_fpath = full_path(gram_dirpath, "encoded_prg");
    parameters.fm_index_fpath = full_path(gram_dirpath, "fm_index");
    parameters.site_mask_fpath = full_path(gram_dirpath, "variant_site_mask");
    parameters.allele_mask_fpath = full_path(gram_dirpath, "allele_mask");
    parameters.sdsl_memory_log_fpath = full_path(gram_dirpath, "sdsl_memory_log");

    parameters.kmer_index_fpath = full_path(gram_dirpath, "kmer_index");
    parameters.kmers_fpath = full_path(gram_dirpath, "kmers");
    parameters.kmer_entry_stats_fpath = full_path(gram_dirpath, "kmer_entry_stats");
    parameters.sa_intervals_fpath = full_path(gram_dirpath, "sa_intervals");
    parameters.paths_fpath = full_path(gram_dirpath, "paths");

    parameters.kmers_size = vm["kmer-size"].as<uint32_t>();
    parameters.max_read_size = vm["max-read-size"].as<uint32_t>();
    return parameters;
}