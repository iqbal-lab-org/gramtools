/* TODO:
 *  *   implement boost logging
*/

#include <iostream>

#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

#include "prg.hpp"
#include "parameters.hpp"
#include "bwt_search.hpp"
#include "map.hpp"
#include "timer_report.hpp"
#include "main.hpp"


int main(int argc, const char *const *argv) {
    auto params = parse_command_line_parameters(argc, argv);
    TimerReport timer_report;

    PRG_Info prg_info;

    std::cout << "Getting FM-index" << std::endl;
    prg_info.fm_index = get_fm_index(true, params.fm_index_fpath, params.prg_integer_alphabet_fpath,
                                           params.linear_prg_fpath, params.fm_index_memory_log_fpath);
    timer_report.record("Construct FM-index");

    std::cout << "Parsing sites_map and allele masks" << std::endl;
    MasksParser masks(params.site_mask_fpath, params.allele_mask_fpath);
    prg_info.sites_mask = masks.sites;
    prg_info.allele_mask = masks.allele;
    prg_info.max_alphabet_num = masks.max_alphabet_num;

    timer_report.record("Parse masks");
    // TODO: should allele_coverage be separated from the masks data structure? No.

    prg_info.dna_rank = calculate_ranks(prg_info.fm_index);
    timer_report.record("Calculating DNA ranks");
    std::cout << "Maximum alphabet number: " << prg_info.max_alphabet_num << std::endl;

    std::cout << "Loading kmer index" << std::endl;
    KmerIndex kmer_index = get_kmer_index(params.prg_kmers_fpath, prg_info);
    timer_report.record("Load kmer index");

    std::cout << "Mapping" << std::endl;
    auto count_mapped = quasimap_reads(kmer_index, masks, params, prg_info.fm_index, prg_info.dna_rank);
    std::cout << "Count mapped: " << count_mapped << std::endl;
    timer_report.record("Mapping");

    std::cout << "Writing allele coverage to file" << std::endl;
    output_allele_coverage(params, masks);
    timer_report.record("Output coverage");

    timer_report.report();
    return 0;
}


Parameters parse_command_line_parameters(int argc, const char *const *argv) {
    namespace po = boost::program_options;
    Parameters params;

    boost::program_options::options_description description("All parameters must be specified");
    description.add_options()
                       ("help", "produce help message")
                       ("prg,marker_porition", po::value<std::string>(&params.linear_prg_fpath),
                        "input file containing linear prg")
                       ("csa,c", po::value<std::string>(&params.fm_index_fpath),
                        "output file where the FM-index is stored")
                       ("input,i", po::value<std::string>(&params.reads_fpath),
                        "reference file (FASTA or FASTQ)")
                       ("ps,s", po::value<std::string>(&params.site_mask_fpath),
                        "input file containing mask over the linear prg that indicates at "
                                "each position whether you are inside a site and if so, which site")
                       ("pa,a", po::value<std::string>(&params.allele_mask_fpath),
                        "input file containing mask over the linear prg that indicates at "
                                "each position whether you are inside a allele and if so, which allele")
                       ("co,v", po::value<std::string>(&params.allele_coverage_fpath),
                        "name of output file where coverages on each allele are printed")
                       ("ro,r", po::value<std::string>(&params.processed_reads_fpath),
                        "name of output file where reads that have been processed are printed")
                       ("po,b", po::value<std::string>(&params.prg_integer_alphabet_fpath),
                        "output filename of binary file containing the prg in integer alphabet")
                       ("log,l", po::value<std::string>(&params.fm_index_memory_log_fpath),
                        "output memory log file for the FM-index")
                       ("kfile,f", po::value<std::string>(&params.prg_kmers_fpath),
                        "input file listing all kmers in PRG")
                       ("ksize,k", po::value<int>(&params.kmers_size),
                        "size of pre-calculated kmers");

    boost::program_options::variables_map vm;
    boost::program_options::store(po::parse_command_line(argc, argv, description), vm);

    try {
        boost::program_options::notify(vm);
    } catch (std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cout << description << std::endl;
        exit(-1);
    }

    if (vm.count("help") != 0) {
        std::cout << description << std::endl;
        exit(0);
    }

    return params;
}
