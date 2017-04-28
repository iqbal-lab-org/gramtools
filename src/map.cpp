/* TODO:
 *  *   implement boost logging
*/

#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

#include <boost/program_options/option.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

#include <seqread.hpp>
#include "bwt_search.h"
#include "parse_masks.h"
#include "precalc_gen.hpp"
#include "map.hpp"


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


uint64_t map_festa(Parameters &params, MasksParser &masks,
                   KmerIdx &kmer_idx, KmerIdx &kmer_idx_rev,
                   KmerSites &kmer_sites, KmersRef &kmers_in_ref, CSA &csa) {
    SeqRead input_festa(params.festa_fpath.c_str());
    std::ofstream reads_fhandle(params.processed_reads_fpath);
    uint64_t count_mapped = 0;
    int count_reads = 0;
    int inc = 0;

    std::vector<uint8_t> readin_integer_seq;
    readin_integer_seq.reserve(200);

    for (auto festa_read: input_festa) {
        if (!(inc++ % 10))
            reads_fhandle << count_reads << std::endl;

        process_festa_sequence(festa_read, readin_integer_seq, params,
                               masks, count_reads, kmer_idx, kmer_idx_rev,
                               kmer_sites, kmers_in_ref, count_mapped, csa);

    }
    reads_fhandle.close();
    return count_mapped;
}


void process_festa_sequence(GenomicRead *festa_read, std::vector<uint8_t> &readin_integer_seq,
                            Parameters &params, MasksParser &masks, int &count_reads,
                            KmerIdx &kmer_idx, KmerIdx &kmer_idx_rev,
                            KmerSites &kmer_sites, KmersRef &kmers_in_ref,
                            uint64_t &count_mapped, CSA &csa) {

    std::cout << festa_read->seq << std::endl;
    bool invalid_base_flag = convert_festa_to_int_seq(festa_read, readin_integer_seq);
    if (invalid_base_flag)
        // TODO: should readin_integer_seq and count_reads be modified here?
        return;

    // TODO: is there a way to avoid making this copy?
    auto kmer_start_it = readin_integer_seq.begin() + readin_integer_seq.size() - params.kmers_size;
    auto kmer_end_it = readin_integer_seq.end();
    std::vector<uint8_t> kmer(kmer_start_it, kmer_end_it);

    bool kmer_within_range = kmer_idx.find(kmer) != kmer_idx.end()
                             && kmer_idx_rev.find(kmer) != kmer_idx_rev.end()
                             && kmer_sites.find(kmer) != kmer_sites.end();
    if (!kmer_within_range){
        count_reads++;
        readin_integer_seq.clear();
        return;
    }

    auto sa_intervals = kmer_idx[kmer];
    auto sa_intervals_rev = kmer_idx_rev[kmer];
    auto sites = kmer_sites[kmer];

    auto it = sa_intervals.begin();
    auto it_rev = sa_intervals_rev.begin();

    // kmers in ref means kmers that do not cross any numbers
    // These are either in non-variable region, or are entirely within alleles
    bool first_del = !(kmers_in_ref.find(kmer) != kmers_in_ref.end());
    bool precalc_done = true;

    bidir_search_bwd(csa, (*it).first, (*it).second,
                     (*it_rev).first, (*it_rev).second,
                     readin_integer_seq.begin(),
                     readin_integer_seq.begin() + readin_integer_seq.size() - params.kmers_size,
                     sa_intervals, sa_intervals_rev,
                     sites, masks.allele, masks.max_alphabet_num, first_del, precalc_done);

    // proxy for mapping is "unique horizontally"
    if (sa_intervals.size() != 1){
        count_reads++;
        readin_integer_seq.clear();
        return;
    }

    if (!first_del)
        // becasue sites has one element with an empty vector
        sites.clear();

    count_mapped++;
    uint64_t no_occ = (*it).second - (*it).first;
    it = sa_intervals.begin();

    for (auto ind = (*it).first; ind < (*it).second; ind++) {
        if (sites.empty()) {
            if (masks.allele[csa[ind]] != 0) {
                masks.allele_coverage[(masks.sites[csa[ind]] - 5) / 2][masks.allele[csa[ind]] - 1]++;
                assert(masks.allele[csa[ind]] == masks.allele[csa[ind] + readin_integer_seq.size() - 1]);
            }
            continue;
        }
        // first_del=true - match in an interval starting with a number, all matches must be just to left of end marker

        // if no_occ>1, two matches both starting at the end marker. If one crossed the start marker,
        // sorina would have split into two SAs and here we are in one.
        // so neither crosses the start marker, both start at the end. Since she only updates sites
        // when you cross the left marker, it should be true that sites.front().back().second.size==0
        if (!(sites.empty()) && (no_occ > 1))
            // vertically non-unique
            assert(sites.front().back().second.size() == 0);

        bool invalid = false;
        for (auto it_s : sites) {
            for (auto site_pair : it_s) {
                auto allele = site_pair.second;
                if (it_s != sites.back() && it_s != sites.front() && allele.empty())
                    invalid = true;
            }
        }

        if (invalid)
            continue;

        for (auto it_s : sites) {
            for (auto site_pair : it_s) {
                auto site = site_pair.first;
                auto allele = site_pair.second;
                if (it_s != sites.back() && it_s != sites.front())
                    assert(allele.size() == 1);
                // mask_a[csa[ind]] can be 0 here if the match is
                // coming from a skipped start_site marker
                if ((allele.empty()) && (masks.allele[csa[ind]] > 0))
                    masks.allele_coverage[(site - 5) / 2][masks.allele[csa[ind]] - 1]++;
                else
                    for (auto al : allele)
                        masks.allele_coverage[(site - 5) / 2][al - 1]++;
            }
        }
    }

    count_reads++;
    readin_integer_seq.clear();
}


bool convert_festa_to_int_seq(GenomicRead *festa_read, std::vector<uint8_t> &readin_integer_seq){
    bool invalid_base_flag = false;
    for (int i = 0; i < strlen(festa_read->seq); i++) {
        if (festa_read->seq[i] == 'A' or festa_read->seq[i] == 'a')
            readin_integer_seq.push_back(1);
        else if (festa_read->seq[i] == 'C' or festa_read->seq[i] == 'c')
            readin_integer_seq.push_back(2);
        else if (festa_read->seq[i] == 'G' or festa_read->seq[i] == 'g')
            readin_integer_seq.push_back(3);
        else if (festa_read->seq[i] == 'T' or festa_read->seq[i] == 't')
            readin_integer_seq.push_back(4);
        else
            // TODO: should there be a break here?
            invalid_base_flag = true;
    }
    return invalid_base_flag;
}


void output_allele_coverage(Parameters &params, MasksParser &masks) {
    std::ofstream allele_coverage_fhandle(params.allele_coverage_fpath);
    for (uint32_t i = 0; i < masks.allele_coverage.size(); i++) {
        for (uint32_t j = 0; j < masks.allele_coverage[i].size(); j++)
            allele_coverage_fhandle << masks.allele_coverage[i][j] << " ";
        allele_coverage_fhandle << std::endl;
    }
}
