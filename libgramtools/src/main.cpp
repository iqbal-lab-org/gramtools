#include <iostream>

#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/variant/variant.hpp>
#include <boost/variant/get.hpp>

#include <sdsl/bit_vectors.hpp>

#include "prg.hpp"
#include "masks.hpp"
#include "timer_report.hpp"
#include "kmer_index.hpp"
#include "coverage_analysis.hpp"
#include "main.hpp"


namespace po = boost::program_options;


int main(int argc, const char *const *argv) {
    Parameters parameters;
    Commands command;
    std::tie(parameters, command) = parse_command_line_parameters(argc, argv);
    switch (command) {
        case Commands::build:
            build(parameters);
            break;
        case Commands::quasimap:
            quasimap(parameters);
            break;
    }
    return 0;
}


void build(const Parameters &parameters) {
    std::cout << "Executing build command" << std::endl;
    auto timer = TimerReport();

    std::cout << "Generating integer encoded PRG" << std::endl;
    timer.start("Encoded PRG");
    auto encoded_prg = generate_encoded_prg(parameters);
    timer.stop();
    std::cout << "Number of charecters in integer encoded linear PRG: "
              << encoded_prg.size()
              << std::endl;

    std::cout << "Generating FM-Index" << std::endl;
    timer.start("Generate FM-Index");
    auto fm_index = generate_fm_index(parameters);
    timer.stop();

    std::cout << "Generating PRG masks" << std::endl;
    timer.start("Generating PRG masks");
    MasksParser masks(parameters.site_mask_fpath);

    auto allele_mask = generate_allele_mask(encoded_prg);
    sdsl::store_to_file(allele_mask, parameters.allele_mask_fpath);

    auto markers_mask = generate_markers_mask(encoded_prg);
    // auto markers_rank = sdsl::rank_support_v<1>(&markers_mask);
    // auto markers_select = sdsl::select_support_mcl<1>(&markers_mask);
    timer.stop();

    PRG_Info prg_info = {
            fm_index,
            encoded_prg,
            masks.sites,
            allele_mask,
            markers_mask,
            masks.max_alphabet_num
    };

    std::cout << "Generating kmer index" << std::endl;
    timer.start("Generate kmer index");
    generate_kmer_index(parameters, prg_info);
    timer.stop();

    timer.report();
}


void quasimap(const Parameters &parameters) {
    std::cout << "Executing quasimap command" << std::endl;
    auto timer = TimerReport();

    std::cout << "Loading data" << std::endl;
    timer.start("Load data");
    const auto prg_info = load_prg_info(parameters);
    const auto kmer_index = load_kmer_index(parameters);
    timer.stop();

    std::cout << "Running quasimap" << std::endl;
    timer.start("Quasimap");
    uint64_t all_reads_count = 0;
    uint64_t skipped_reads_count = 0;
    uint64_t mapped_reads_count = 0;

    std::tie(all_reads_count,
             skipped_reads_count,
             mapped_reads_count) = quasimap_reads(parameters,
                                                  kmer_index,
                                                  prg_info);
    std::cout << "Count all reads: " << all_reads_count << std::endl;
    std::cout << "Count skipped reads: " << skipped_reads_count << std::endl;
    std::cout << "Count mapped reads: " << mapped_reads_count << std::endl;
    timer.stop();

    timer.report();
}


Parameters parse_build_parameters(po::variables_map &vm, const po::parsed_options &parsed) {
    po::options_description build_description("build options");
    build_description.add_options()
                             ("prg", po::value<std::string>(),
                              "file containing a linear PRG")
                             ("encoded-prg", po::value<std::string>(),
                              "output file containing the integer encoded linear PRG")
                             ("fm-index", po::value<std::string>(),
                              "file containing the SDSL FM-Index")
                             ("variant-site-mask", po::value<std::string>(),
                              "vaiant site boundary marker mask over the linear PRG")
                             ("allele-mask", po::value<std::string>(),
                              "allele ID mask over the linear PRG")
                             ("memory-log", po::value<std::string>(),
                              "SDSL memory log of FM-index construction")
                             ("kmers-prefix-diffs", po::value<std::string>(),
                              "file continaing kmer prefix diffs used in consturcting the kmer index")
                             ("kmer-index", po::value<std::string>(),
                              "output destination of the kmer index")
                             ("kmer-size", po::value<uint32_t>(),
                              "kmer size used in constructing the kmer index");

    std::vector<std::string> opts = po::collect_unrecognized(parsed.options,
                                                             po::include_positional);
    opts.erase(opts.begin());

    po::store(po::command_line_parser(opts).options(build_description).run(), vm);

    Parameters parameters;
    parameters.linear_prg_fpath = vm["prg"].as<std::string>();
    parameters.encoded_prg_fpath = vm["encoded-prg"].as<std::string>();
    parameters.fm_index_fpath = vm["fm-index"].as<std::string>();
    parameters.site_mask_fpath = vm["variant-site-mask"].as<std::string>();
    parameters.allele_mask_fpath = vm["allele-mask"].as<std::string>();
    parameters.sdsl_memory_log_fpath = vm["memory-log"].as<std::string>();
    parameters.kmer_suffix_diffs_fpath = vm["kmers-prefix-diffs"].as<std::string>();
    parameters.kmer_index_fpath = vm["kmer-index"].as<std::string>();
    parameters.kmers_size = vm["kmer-size"].as<uint32_t>();
    return parameters;
}


Parameters parse_quasimap_parameters(po::variables_map &vm,
                                     const po::parsed_options &parsed) {
    po::options_description quasimap_description("quasimap options");
    quasimap_description.add_options()
                                ("prg", po::value<std::string>(),
                                 "file containing a linear PRG")
                                ("encoded-prg", po::value<std::string>(),
                                 "output file containing the integer encoded linear PRG")
                                ("fm-index", po::value<std::string>(),
                                 "file containing the SDSL FM-Index")
                                ("variant-site-mask", po::value<std::string>(),
                                 "vaiant site boundary marker mask over the linear PRG")
                                ("allele-mask", po::value<std::string>(),
                                 "allele ID mask over the linear PRG")
                                ("memory-log", po::value<std::string>(),
                                 "SDSL memory log of FM-index construction")
                                ("kmers-prefix-diffs", po::value<std::string>(),
                                 "file continaing kmer prefix diffs used in consturcting the kmer index")
                                ("kmer-index", po::value<std::string>(),
                                 "output destination of the kmer index")
                                ("kmer-size", po::value<uint32_t>(),
                                 "kmer size used in constructing the kmer index")

                                ("reads", po::value<std::string>(),
                                 "file contining reads (FASTA or FASTQ)")
                                ("allele-coverages", po::value<std::string>(),
                                 "output file of read coverages over each allele and variant site")
                                ("reads-progress", po::value<std::string>(),
                                 "output file containing progress counts of reads processed");

    std::vector<std::string> opts = po::collect_unrecognized(parsed.options,
                                                             po::include_positional);
    opts.erase(opts.begin());
    po::store(po::command_line_parser(opts).options(quasimap_description).run(), vm);

    Parameters parameters;
    parameters.linear_prg_fpath = vm["prg"].as<std::string>();
    parameters.encoded_prg_fpath = vm["encoded-prg"].as<std::string>();
    parameters.fm_index_fpath = vm["fm-index"].as<std::string>();
    parameters.site_mask_fpath = vm["variant-site-mask"].as<std::string>();
    parameters.allele_mask_fpath = vm["allele-mask"].as<std::string>();
    parameters.sdsl_memory_log_fpath = vm["memory-log"].as<std::string>();
    parameters.kmer_suffix_diffs_fpath = vm["kmers-prefix-diffs"].as<std::string>();
    parameters.kmer_index_fpath = vm["kmer-index"].as<std::string>();
    parameters.kmers_size = vm["kmer-size"].as<uint32_t>();

    parameters.reads_fpath = vm["reads"].as<std::string>();
    parameters.allele_coverage_fpath = vm["allele-coverages"].as<std::string>();
    parameters.reads_progress_fpath = vm["reads-progress"].as<std::string>();

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
        auto parameters = parse_quasimap_parameters(vm, parsed);
        return std::make_pair(parameters, Commands::quasimap);
    }

    // unrecognised command
    throw po::invalid_option_value(cmd);
}
