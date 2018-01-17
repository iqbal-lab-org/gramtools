#include <iostream>

#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/variant/variant.hpp>
#include <boost/variant/get.hpp>
#include <boost/filesystem.hpp>

#include <sdsl/bit_vectors.hpp>

#include "prg/prg.hpp"
#include "prg/masks.hpp"
#include "timer_report.hpp"
#include "kmer_index/kmer_index.hpp"
#include "quasimap/quasimap.hpp"
#include "main.hpp"


namespace po = boost::program_options;
namespace fs = boost::filesystem;


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

    PRG_Info prg_info;

    std::cout << "Generating integer encoded PRG" << std::endl;
    timer.start("Encoded PRG");
    prg_info.encoded_prg = generate_encoded_prg(parameters);
    timer.stop();
    std::cout << "Number of charecters in integer encoded linear PRG: "
              << prg_info.encoded_prg.size()
              << std::endl;

    std::cout << "Generating FM-Index" << std::endl;
    timer.start("Generate FM-Index");
    prg_info.fm_index = generate_fm_index(parameters);
    timer.stop();

    std::cout << "Generating PRG masks" << std::endl;
    timer.start("Generating PRG masks");
    generate_dna_bwt_masks(prg_info.fm_index, parameters);

    MasksParser masks(parameters.site_mask_fpath);
    prg_info.sites_mask = masks.sites;
    prg_info.max_alphabet_num = masks.max_alphabet_num;

    prg_info.allele_mask = generate_allele_mask(prg_info.encoded_prg);
    sdsl::store_to_file(prg_info.allele_mask, parameters.allele_mask_fpath);

    prg_info.prg_markers_mask = generate_prg_markers_mask(prg_info.encoded_prg);
    prg_info.prg_markers_rank = sdsl::rank_support_v<1>(&prg_info.prg_markers_mask);
    prg_info.prg_markers_select = sdsl::select_support_mcl<1>(&prg_info.prg_markers_mask);

    prg_info.bwt_markers_mask = generate_bwt_markers_mask(prg_info.fm_index);
    prg_info.bwt_markers_rank = sdsl::rank_support_v<1>(&prg_info.bwt_markers_mask);
    prg_info.bwt_markers_select = sdsl::select_support_mcl<1>(&prg_info.bwt_markers_mask);
    prg_info.markers_mask_count_set_bits =
            prg_info.bwt_markers_rank(prg_info.bwt_markers_mask.size());

    prg_info.dna_bwt_masks = load_dna_bwt_masks(prg_info.fm_index, parameters);
    prg_info.rank_bwt_a = sdsl::rank_support_v<1>(&prg_info.dna_bwt_masks.mask_a);
    prg_info.rank_bwt_c = sdsl::rank_support_v<1>(&prg_info.dna_bwt_masks.mask_c);
    prg_info.rank_bwt_g = sdsl::rank_support_v<1>(&prg_info.dna_bwt_masks.mask_g);
    prg_info.rank_bwt_t = sdsl::rank_support_v<1>(&prg_info.dna_bwt_masks.mask_t);
    timer.stop();

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
    auto quasimap_stats = quasimap_reads(parameters, kmer_index, prg_info);
    std::cout << "Count all reads: " << quasimap_stats.all_reads_count << std::endl;
    std::cout << "Count skipped reads: " << quasimap_stats.skipped_reads_count << std::endl;
    std::cout << "Count mapped reads: " << quasimap_stats.mapped_reads_count << std::endl;
    timer.stop();

    timer.report();
}


std::string full_path(const std::string &gram_dirpath,
                      const std::string &file_name) {
    fs::path dir(gram_dirpath);
    fs::path file(file_name);
    fs::path full_path = dir / file;
    return full_path.string();
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


Parameters parse_quasimap_parameters(po::variables_map &vm,
                                     const po::parsed_options &parsed) {
    po::options_description quasimap_description("quasimap options");
    quasimap_description.add_options()
                                ("gram", po::value<std::string>(),
                                 "gramtools directory")
                                ("reads", po::value<std::string>(),
                                 "file contining reads (FASTA or FASTQ)")
                                ("kmer-size", po::value<uint32_t>(),
                                 "kmer size used in constructing the kmer index")
                                ("run-directory", po::value<std::string>(),
                                 "a directory which contains all quasimap output files");

    std::vector<std::string> opts = po::collect_unrecognized(parsed.options,
                                                             po::include_positional);
    opts.erase(opts.begin());
    po::store(po::command_line_parser(opts).options(quasimap_description).run(), vm);

    std::string gram_dirpath = vm["gram"].as<std::string>();

    Parameters parameters;
    parameters.gram_dirpath = gram_dirpath;
    parameters.linear_prg_fpath = full_path(gram_dirpath, "prg");
    parameters.encoded_prg_fpath = full_path(gram_dirpath, "encoded_prg");
    parameters.fm_index_fpath = full_path(gram_dirpath, "fm_index");
    parameters.site_mask_fpath = full_path(gram_dirpath, "variant_site_mask");
    parameters.allele_mask_fpath = full_path(gram_dirpath, "allele_mask");
    parameters.kmer_index_fpath = full_path(gram_dirpath, "kmer_index");
    
    parameters.kmers_size = vm["kmer-size"].as<uint32_t>();
    parameters.reads_fpath = vm["reads"].as<std::string>();

    std::string run_dirpath = vm["run-directory"].as<std::string>();

    parameters.allele_coverage_fpath = full_path(run_dirpath, "allele_sum_coverage");
    parameters.reads_progress_fpath = full_path(run_dirpath, "reads_progress");
    parameters.sdsl_memory_log_fpath = full_path(run_dirpath, "sdsl_memory_log");

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
