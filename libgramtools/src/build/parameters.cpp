#include <iostream>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

#include "common/parameters.hpp"

#include "build/parameters.hpp"

using namespace gram;
using namespace gram::commands::build;


BuildParams commands::build::parse_parameters(po::variables_map &vm, const po::parsed_options &parsed) {

    std::string gram_dirpath, fasta_ref;
    uint32_t kmer_size;
    uint32_t max_read_size;

    po::options_description build_description("build options");
    build_description.add_options()
            ("gram_dir", po::value<std::string>(&gram_dirpath)->required(),
             "gramtools directory")
            ("ref", po::value<std::string>(&fasta_ref)->required(),
             "fasta reference file")
            ("kmer_size", po::value<uint32_t>(&kmer_size)->required(),
             "kmer size used in constructing the kmer index")
            ("max_threads", po::value<uint32_t>()->default_value(1),
             "maximum number of threads used")
            ("all_kmers", po::bool_switch()->default_value(false),
             "generate all kmers of given size (as opposed to inspecting PRG for min set)")
            ("max_read_size", po::value<uint32_t>(&max_read_size)->default_value(0),
                    "read maximum size for the set of reads used when quasimaping");

    std::vector<std::string> opts = po::collect_unrecognized(parsed.options,
                                                             po::include_positional);
    if (opts.size() > 0) opts.erase(opts.begin()); // Takes out the command itself
    try {
        po::store(po::command_line_parser(opts).options(build_description).run(), vm);
        po::notify(vm);
    }
    catch(const std::exception &e){
        std::cerr << e.what() << std::endl;
        std::cerr << build_description << std::endl;
        exit(1);
    }

    BuildParams parameters = {};
    fill_common_parameters(parameters, gram_dirpath);

    parameters.sdsl_memory_log_fpath = full_path(gram_dirpath, "sdsl_memory_log");
    parameters.kmers_size = kmer_size;
    parameters.fasta_ref = fasta_ref;

    parameters.all_kmers_flag = vm["all_kmers"].as<bool>();
    parameters.max_read_size = vm["max_read_size"].as<uint32_t>();
    parameters.maximum_threads = vm["max_threads"].as<uint32_t>();

    if (! parameters.all_kmers_flag and parameters.max_read_size == 0)
            throw std::invalid_argument(
                    "--max_read_size must be > 0 when --all_kmers flag is not used"
                    );

    return parameters;
}