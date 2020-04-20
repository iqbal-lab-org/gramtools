#include <iostream>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

#include "common/parameters.hpp"

#include "build/parameters.hpp"

using namespace gram;
using namespace gram::commands::build;


BuildParams commands::build::parse_parameters(po::variables_map &vm, const po::parsed_options &parsed) {

    std::string gram_dirpath;
    uint32_t kmer_size;
    uint32_t max_read_size;

    po::options_description build_description("build options");
    build_description.add_options()
                             ("gram", po::value<std::string>(&gram_dirpath)->required(),
                              "gramtools directory")
                             ("kmer-size", po::value<uint32_t>(&kmer_size)->required(),
                              "kmer size used in constructing the kmer index")
                             ("max-read-size", po::value<uint32_t>(&max_read_size)->required(),
                              "read maximum size for the set of reads used when quasimaping")
                             ("max-threads", po::value<uint32_t>()->default_value(1),
                              "maximum number of threads used")
                             ("all-kmers", po::bool_switch()->default_value(false),
                              "generate all kmers of given size (as opposed to inspecting PRG for min set)");

    std::vector<std::string> opts = po::collect_unrecognized(parsed.options,
                                                             po::include_positional);
    if (opts.size() > 0) opts.erase(opts.begin()); // Takes out the command itself
    try {
        po::store(po::command_line_parser(opts).options(build_description).run(), vm);
        po::notify(vm);
    }
    catch(const std::exception &e){
        std::cout << e.what() << std::endl;
        std::cout << build_description << std::endl;
        exit(1);
    }

    BuildParams parameters = {};
    fill_common_parameters(parameters, gram_dirpath);

    parameters.sdsl_memory_log_fpath = full_path(gram_dirpath, "sdsl_memory_log");
    parameters.kmers_size = kmer_size;
    parameters.max_read_size = max_read_size;

    parameters.all_kmers_flag = vm["all-kmers"].as<bool>();
    parameters.maximum_threads = vm["max-threads"].as<uint32_t>();

    return parameters;
}