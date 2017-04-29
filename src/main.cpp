/* TODO:
 *  *   implement boost logging
*/

#include <fstream>
#include <iostream>
#include <vector>

#include <boost/program_options/option.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

#include <seqread.hpp>
#include "parameters.hpp"
#include "bwt_search.h"
#include "masks.hpp"
#include "kmers.hpp"
#include "map.hpp"
#include "main.hpp"


namespace po = boost::program_options;


int main(int argc, const char *const *argv) {
    auto params = parse_command_line_parameters(argc, argv);

    std::cout << "Constructing CSA" << std::endl;
    CSA csa = csa_constr(params.prg_fpath, params.prg_integer_alphabet_fpath,
                         params.csa_memory_log_fpath, params.csa_fpath, true, true);

    std::cout << "Parsing sites and allele masks" << std::endl;
    MasksParser masks(params.site_mask_fpath, params.allele_mask_fpath);
    // TODO: should allele_coverage be separated from the masks data structure?

    std::cout << "Pre-calculating K-mers" << std::endl;
    KmerIdx kmer_idx, kmer_idx_rev;
    KmerSites kmer_sites;
    KmersRef kmers_in_ref;
    get_precalc_kmers(csa, kmer_idx, kmer_idx_rev,
                      kmer_sites, kmers_in_ref, masks.allele,
                      params.prg_kmers_fpath, masks.max_alphabet_num, params.kmers_size);

    std::cout << "Mapping" << std::endl;
    uint64_t count_mapped = map_festa(params, masks, kmer_idx, kmer_idx_rev,
                                      kmer_sites, kmers_in_ref, csa);
    std::cout << "Count mapped: " << count_mapped << std::endl;

    std::cout << "Writing allele coverage to file" << std::endl;
    output_allele_coverage(params, masks);
    return 0;
}


Parameters parse_command_line_parameters(int argc, const char *const *argv) {
    Parameters params;
    po::options_description description("All parameters must be specified");
    description.add_options()
            ("help", "produce help message")
            ("prg,p", po::value<std::string>(&params.prg_fpath)->required(),
             "input file containing linear prg")
            ("csa,c", po::value<std::string>(&params.csa_fpath)->required(),
             "output file where CSA is stored")
            ("input,i", po::value<std::string>(&params.festa_fpath)->required(),
             "input FASTA/FASTQ file to be mapped")
            ("ps,s", po::value<std::string>(&params.site_mask_fpath)->required(),
             "input file containing mask over the linear prg that indicates at "
                     "each position whether you are inside a site and if so, which site")
            ("pa,c", po::value<std::string>(&params.allele_mask_fpath)->required(),
             "input file containing mask over the linear prg that indicates at "
                     "each position whether you are inside a allele and if so, which allele")
            ("co,v", po::value<std::string>(&params.allele_coverage_fpath)->required(),
             "name of output file where coverages on each allele are printed")
            ("ro,r", po::value<std::string>(&params.processed_reads_fpath)->required(),
             "name of output file where reads that have been processed are printed")
            ("po,b", po::value<std::string>(&params.prg_integer_alphabet_fpath)->required(),
             "output filename of binary file containing the prg in integer alphabet")
            ("log,l", po::value<std::string>(&params.csa_memory_log_fpath)->required(),
             "output memory log file for CSA")
            ("kfile,f", po::value<std::string>(&params.prg_kmers_fpath)->required(),
             "input file listing all kmers in PRG")
            ("ksize,k", po::value<int>(&params.kmers_size)->required(),
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
