#include "simulate/parameters.hpp"
#include <iostream>

using namespace gram;

SimulateParams commands::simulate::parse_parameters(po::variables_map &vm,
                                                    const po::parsed_options &parsed){

    SimulateParams parameters;
    std::string output_dir_fpath;

    po::options_description simulate_description("simulate options");
    simulate_description.add_options()
            ("prg", po::value<std::string>(&parameters.encoded_prg_fpath)->required(),
             "path to encoded prg")
            ("n", po::value<uint64_t>(&parameters.max_num_paths)->required(),
             "Max num paths through prg to simulate")
            ("sample_id", po::value<std::string>(&parameters.sample_id)->required(),
             "Prefixes output filenames and sample IDs in output files")
            ("o", po::value<std::string>(&output_dir_fpath)->required(),
             "directory containing outputs")
            ("i", po::value<std::string>(&parameters.input_sequences_fpath),
                    "input sequences to induce genotypes on");

    std::vector<std::string> opts = po::collect_unrecognized(parsed.options,
                                                             po::include_positional);
    if (opts.size() > 0) opts.erase(opts.begin()); // Takes out the command itself
    try {
        po::store(po::command_line_parser(opts).options(simulate_description).run(), vm);
        po::notify(vm);
    }
    catch(const std::exception &e){
        std::cout << e.what() << std::endl;
        std::cout << simulate_description << std::endl;
        exit(1);
    }

    parameters.json_out_fpath = full_path(output_dir_fpath,
            parameters.sample_id + std::string(".json"));
    parameters.fasta_out_fpath = full_path(output_dir_fpath,
            parameters.sample_id + std::string(".fasta"));
    return parameters;
}
