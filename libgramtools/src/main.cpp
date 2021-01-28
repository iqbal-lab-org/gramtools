#include <iostream>
#include <sdsl/bit_vectors.hpp>

#include "build/build.hpp"
#include "build/parameters.hpp"
#include "common/timer_report.hpp"
#include "genotype/genotype.hpp"
#include "genotype/parameters.hpp"
#include "prg/prg_info.hpp"
#include "simulate/parameters.hpp"
#include "simulate/simulate.hpp"

using namespace gram;

namespace gram {
enum class Command { build, genotype, simulate };

struct top_level_params {
  po::variables_map vm;
  po::parsed_options parsed;
  Command command;
};

top_level_params parse_command_line_parameters(int argc,
                                               const char *const *argv);
}  // namespace gram

int main(int argc, const char *const *argv) {
  auto command_params = parse_command_line_parameters(argc, argv);
  auto debug = command_params.vm["debug"].as<bool>();

  if (command_params.command == Command::build) {
    BuildParams build_params = commands::build::parse_parameters(
        command_params.vm, command_params.parsed);
    commands::build::run(build_params);
  } else if (command_params.command == Command::genotype) {
    GenotypeParams geno_params = commands::genotype::parse_parameters(
        command_params.vm, command_params.parsed);
    commands::genotype::run(geno_params, debug);
  }

  else if (command_params.command == Command::simulate) {
    SimulateParams simu_params = commands::simulate::parse_parameters(
        command_params.vm, command_params.parsed);
    commands::simulate::run(simu_params);
  }

  return 0;
}

top_level_params gram::parse_command_line_parameters(int argc,
                                                     const char *const *argv) {
  po::options_description global("Gramtools! Global options");
  global.add_options()("command", po::value<std::string>(),
                       "command to execute: {build, genotype, simulate}")(
      "subargs", po::value<std::vector<std::string> >(),
      "arguments to command")("help", "Produce this help message")(
      "debug", po::bool_switch()->default_value(false), "Turn on debug output");

  // We register all positional arguments after the command as 'subargs'
  // which get forwarded to the command specific parser.
  po::positional_options_description pos;
  pos.add("command", 1).add("subargs", -1);

  po::variables_map vm;

  po::parsed_options parsed =
      po::command_line_parser(argc, argv)
          .options(global)
          .positional(pos)
          .allow_unregistered()  // Allows passing of subargs to command
                                 // specific parser.
          .run();

  po::store(parsed, vm);

  // Note: could use po::notify() and make 'command' required too- see command
  // specific parsers.
  if (vm.count("help") || !vm.count("command")) {
    std::cout << global << std::endl;
    exit(0);
  }

  std::string cmd_string = vm["command"].as<std::string>();

  Command cmd;
  if (cmd_string == "build")
    cmd = Command::build;
  else if (cmd_string == "genotype")
    cmd = Command::genotype;
  else if (cmd_string == "simulate")
    cmd = Command::simulate;
  else {
    std::cout << "Unrecognised command: " << cmd_string << std::endl;
    std::cout << global << std::endl;
    exit(1);
  }

  return top_level_params{vm, parsed, cmd};
}
