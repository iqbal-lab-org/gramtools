#include <iostream>

#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/variant/variant.hpp>
#include <boost/variant/get.hpp>
#include <boost/filesystem.hpp>

#include <sdsl/bit_vectors.hpp>

#include "prg/prg.hpp"
#include "common/timer_report.hpp"

#include "build/build.hpp"
#include "quasimap/quasimap.hpp"
#include "quasimap/parameters.hpp"

#include "main.hpp"


namespace po = boost::program_options;
namespace fs = boost::filesystem;


int main(int argc, const char *const *argv) {
    Parameters parameters;
    Commands command;
    std::tie(parameters, command) = parse_command_line_parameters(argc, argv);
    switch (command) {
        case Commands::build:
            commands::build::run(parameters);
            break;
        case Commands::quasimap:
            commands::quasimap::run(parameters);
            break;
    }
    return 0;
}


Parameters parse_build_parameters(po::variables_map &vm, const po::parsed_options &parsed) {
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
    parameters.kmer_index_fpath = full_path(gram_dirpath, "kmer_index");
    parameters.sdsl_memory_log_fpath = full_path(gram_dirpath, "sdsl_memory_log");

    parameters.kmers_size = vm["kmer-size"].as<uint32_t>();
    parameters.max_read_size = vm["max-read-size"].as<uint32_t>();
    return parameters;
}


std::pair<Parameters, Commands> parse_command_line_parameters(int argc, const char *const *argv) {
    po::options_description global("Global options");
    global.add_options()
                  ("debug", "Turn on debug output")
                  ("command", po::value<std::string>(), "command to execute")
                  ("subargs", po::value<std::vector<std::string> >(), "Arguments for command");

    po::positional_options_description pos;
    pos.add("command", 1)
       .add("subargs", -1);

    po::variables_map vm;

    po::parsed_options parsed =
            po::command_line_parser(argc, argv)
                    .options(global)
                    .positional(pos)
                    .allow_unregistered()
                    .run();

    po::store(parsed, vm);
    std::string cmd = vm["command"].as<std::string>();

    if (cmd == "build") {
        auto parameters = parse_build_parameters(vm, parsed);
        return std::make_pair(parameters, Commands::build);
    } else if (cmd == "quasimap") {
        auto parameters = commands::quasimap::parse_parameters(vm, parsed);
        return std::make_pair(parameters, Commands::quasimap);
    }

    // unrecognised command
    throw po::invalid_option_value(cmd);
}
