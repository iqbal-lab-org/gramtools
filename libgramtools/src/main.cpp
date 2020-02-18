#include <iostream>


#include <sdsl/bit_vectors.hpp>
#include <omp.h>

#include "prg/prg_info.hpp"
#include "common/timer_report.hpp"

#include "build/build.hpp"
#include "build/parameters.hpp"

#include "genotype/genotype.hpp"
#include "genotype/parameters.hpp"


using namespace gram;

namespace gram {
    enum class Command {
        build,
        genotype
    };

    struct top_level_params{
        po::variables_map vm;
        po::parsed_options parsed;
        Command command;
    };

   top_level_params parse_command_line_parameters(int argc, const char *const *argv);
}

int main(int argc, const char *const *argv) {
    auto command_params = parse_command_line_parameters(argc, argv);

    CommonParameters common_params;
    BuildParams build_params;
    GenotypeParams geno_params;

    switch(command_params.command){
        case Command::build:
            build_params =
                    commands::build::parse_parameters(command_params.vm, command_params.parsed);
            common_params = build_params;
            commands::build::run(build_params);
            break;

        case Command::genotype:
            geno_params =
                    commands::genotype::parse_parameters(command_params.vm, command_params.parsed);
            common_params = geno_params;
            commands::genotype::run(geno_params);
            break;
    }

    std::cout << "maximum thread count: " << common_params.maximum_threads << std::endl;
    omp_set_num_threads(common_params.maximum_threads);
    return 0;
}


top_level_params gram::parse_command_line_parameters(int argc, const char *const *argv) {
    po::options_description global("Gramtools! Global options");
    global.add_options()
                  ("command", po::value<std::string>(), "command to execute: {build, genotype}")
                  ("subargs", po::value<std::vector<std::string> >(), "arguments to command")
                  ("help", "Produce this help message")
                ("debug", "Turn on debug output");

    // We register all positional arguments after the command as 'subargs'
    // which get forwarded to the command specific parser.
    po::positional_options_description pos;
    pos.add("command", 1)
       .add("subargs", -1);

    po::variables_map vm;

    po::parsed_options parsed = po::command_line_parser(argc, argv)
                    .options(global)
                    .positional(pos)
                    .allow_unregistered() // Allows passing of subargs to command specific parser.
                    .run();

    po::store(parsed, vm);

    //Note: could use po::notify() and make 'command' required too- see command specific parsers.
    if (vm.count("help") || !vm.count("command")) {
        std::cout << global << std::endl;
        exit(1);
    }

    std::string cmd_string = vm["command"].as<std::string>();

    Command cmd;
    if (cmd_string == "build") cmd = Command::build;
    else if (cmd_string == "genotype") cmd = Command::genotype;
    else {
        std::cout << "Unrecognised command: " << cmd_string << std::endl;
        std::cout << global << std::endl;
        exit(1);
    }

    return top_level_params{vm, parsed, cmd};
}
