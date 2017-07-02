/* TODO:
 *  *   implement boost logging
*/

#include <iostream>

#include <boost/timer/timer.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

#include "git_version/git_version.hpp"

#include "parameters.hpp"
#include "bwt_search.h"
#include "masks.hpp"
#include "kmers.hpp"
#include "map.hpp"
#include "process_prg.hpp"
#include "main.hpp"


int main(int argc, const char *const *argv) {
    auto params = parse_command_line_parameters(argc, argv);
    TimerReport timer_report;

    std::cout << "Constructing FM-index" << std::endl;
    FM_Index fm_index = construct_fm_index(params.prg_fpath,
                                           params.prg_integer_alphabet_fpath,
                                           params.fm_index_memory_log_fpath,
                                           params.fm_index_fpath, true);
    timer_report.record("Construct FM-index");

    std::cout << "Parsing sites and allele masks" << std::endl;
    MasksParser masks(params.site_mask_fpath, params.allele_mask_fpath);
    timer_report.record("Parse masks");
    // TODO: should allele_coverage be separated from the masks data structure?
    std::cout << "Maximum alphabet number: " << masks.max_alphabet_num << std::endl;
    
    std::cout << "Pre-calculating kmers" << std::endl;
    KmersData kmers = get_kmers(fm_index, masks.allele, params.prg_kmers_fpath,
                                masks.max_alphabet_num, params.kmers_size);
    timer_report.record("Pre-calc kmers");

    std::cout << "Mapping" << std::endl;
    int no_mapped, count_attempt_mapped;
    std::tie (no_mapped, count_attempt_mapped) = map_festa(params, masks, kmers, fm_index);
    std::cout << "Count mapped: " << no_mapped << " out of "<< count_attempt_mapped << std::endl;
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

    po::options_description description("All parameters must be specified");
    description.add_options()
            ("help", "produce help message")
            ("prg,marker_porition", po::value<std::string>(&params.prg_fpath),
             "input file containing linear prg")
            ("csa,c", po::value<std::string>(&params.fm_index_fpath),
             "output file where CSA is stored")
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
             "output memory log file for CSA")
            ("kfile,f", po::value<std::string>(&params.prg_kmers_fpath),
             "input file listing all kmers in PRG")
            ("ksize,k", po::value<int>(&params.kmers_size),
             "size of pre-calculated kmers");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, description), vm);

    try {
        po::notify(vm);
    } catch (std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cout << description << std::endl;
        exit(-1);
    }

    if (vm.count("help")) {
        std::cout << description << std::endl;
        exit(0);
    }

    return params;
}


void TimerReport::record(std::string note){
    boost::timer::cpu_times times = timer.elapsed();
    double elapsed_time = (times.user + times.system) * 1e-9;
    Entry entry = std::make_pair(note, elapsed_time);
    logger.push_back(entry);
}


void TimerReport::report() const{
    std::cout << "\nTimer report:" << std::endl;
    cout_row(" ", "seconds");

    for (const auto &entry: TimerReport::logger){
        auto &note = std::get<0>(entry);
        auto &elapsed_time = std::get<1>(entry);
        cout_row(note, elapsed_time);
    }
}


template <typename TypeCol1, typename TypeCol2>
void TimerReport::cout_row(TypeCol1 col1, TypeCol2 col2) const {
    std::cout << std::setw(20) << std::right << col1
              << std::setw(10) << std::right << col2
              << std::endl;
}
