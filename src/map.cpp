#include <getopt.h>
#include <time.h>

#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

#include <seqread.hpp>
#include "bwt_search.h"
#include "precalc_gen.hpp"


void timestamp() {
    time_t ltime;
    ltime = time(NULL);
    printf("\n-----\n%s", asctime(localtime(&ltime)));
    fflush(stdout);
}


const std::string usage_statment =
        "\ngramtools usage:\n"
                "All paramaters must be specified.\n\n"
                "--prg    -p   input file containing linear prg\n"
                "--csa    -c   output file where CSA is stored\n"
                "--input  -i   input FASTA/FASTQ file to be mapped\n"
                "--ps     -s   input file containing mask over the\n\t\tlinear prg that indicates at each\n\t\tposition whether you are inside a\n\t\tsite and if so, which site\n"
                "--pa     -a   input file containing mask over the\n\t\tlinear prg that indicates at each\n\t\tposition whether you are inside a\n\t\tsite and if so, which allele\n"
                "--co     -v   name of output file where coverages on each allele are printed\n"
                "--ro     -r   name of output file where reads that have been processed are printed\n"
                "--po     -b   output filename of binary file containing the prg in integer alphabet\n"
                "--log    -l   Output memory log file for CSA\n"
                "--ksize  -k   size of precalculated kmers\n"
                "--kfile  -f   input file listing all kmers in PRG\n";


/**
	argv[1] -  file containing linear prg
	argv[2] -  file where CSA is stored
	argv[3] -  file containing reads to be mapped (one read per line)
	argv[4] -  file containing mask over the linear prg that indicates at each position whether you are inside a site and if so, which site
	argv[5] -  file containing mask over the linear prg that indicates at each position whether you are inside a site and if so, which allele
	argv[6] -  name of output file where coverages on each allele are printed
	argv[7] -  name of output file where reads that have been processed are printed
	argv[8] -  name of binary file where the prg in integer alphabet is written
	argv[9] -  memory log file for CSA
	argv[10] - size of pre-calculated kmers
	argv[11] - kmer file
*/
int main(int argc, char *argv[]) {
    std::string linear_prg_fname, csa_fname, fasta_fname,
            site_mask_fname, allele_mask_fname,
            allele_coverage_fname, reads_fname, out_prg_fname,
            memory_log_fname, input_kmers_size, kmer_fname;

    std::vector<std::string *> pars = {
            &linear_prg_fname, &csa_fname, &fasta_fname,
            &site_mask_fname, &allele_mask_fname,
            &allele_coverage_fname, &reads_fname, &out_prg_fname,
            &memory_log_fname, &input_kmers_size, &kmer_fname
    };

    while (1) {
        static struct option long_options[] =
                {
                        {"prg",   required_argument, 0, 'p'},
                        {"csa",   required_argument, 0, 'c'},
                        {"input", required_argument, 0, 'i'},
                        {"ps",    required_argument, 0, 's'},
                        {"pa",    required_argument, 0, 'a'},
                        {"co",    required_argument, 0, 'v'},
                        {"ro",    required_argument, 0, 'r'},
                        {"po",    required_argument, 0, 'b'},
                        {"log",   required_argument, 0, 'l'},
                        {"ksize", required_argument, 0, 'k'},
                        {"kfile", required_argument, 0, 'f'},
                        {0, 0,                       0, 0}
                };

        // getopt_long stores the option index here.
        int option_index = 0;
        int c = getopt_long(argc, argv, "p:c:i:s:a:b:v:r:l:k:f:", long_options, &option_index);
        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c) {
            case 'p':
            case 'c':
            case 'i':
            case 's':
            case 'a':
            case 'b':
            case 'v':
            case 'r':
            case 'l':
            case 'k':
            case 'f':
                for (unsigned int i = 0; i < pars.size(); i++) {
                    if (long_options[i].val == c)
                        *pars[i] = optarg;
                }
                break;

            case '?':
                // Error message already printed by getopt_long.
                std::cout << "Error parsing arguments\n" << usage_statment;
                break;

            default:
                std::cout << "Error parsing arguments\n" << usage_statment;
                abort();
        }
    }

    for (auto i : pars) {
        if (*i == "") {
            std::cout << "All paramaters not specified.\n" << usage_statment;
            exit(-1);
        }
    }

    timestamp();
    cout << "Start CSA construction" << endl;
    auto csa = csa_constr(linear_prg_fname, out_prg_fname,
                          memory_log_fname, csa_fname, true, true);
    timestamp();
    cout << "End CSA construction" << endl;

    std::vector<uint64_t> mask_sites;
    std::vector<int> mask_allele;
    std::vector<std::vector<int> > covgs;
    uint64_t maxx = parse_masks(mask_sites, mask_allele, site_mask_fname, allele_mask_fname, covgs);

    std::vector<std::vector<string> > site_reads(covgs.size(), std::vector<string>(1));

    int kmers_size = std::atoi(input_kmers_size.c_str());
    sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t, uint64_t>>> kmer_idx, kmer_idx_rev;
    sequence_map<std::vector<uint8_t>, std::list<std::vector<std::pair<uint32_t, std::vector<int>>>>> kmer_sites;
    sequence_set<std::vector<uint8_t>> kmers_in_ref;
    get_precalc_kmers(csa, kmer_idx, kmer_idx_rev, kmer_sites, kmers_in_ref, mask_allele, kmer_fname, maxx, kmers_size);

    cout << "Start mapping" << endl;
    timestamp();

    std::list<std::pair<uint64_t, uint64_t>> sa_intervals, sa_intervals_rev;
    std::list<std::pair<uint64_t, uint64_t>>::iterator it, it_rev;
    std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;
    std::list<std::vector<std::pair<uint32_t, std::vector<int>>>>::iterator it_s;

    std::vector<uint8_t> readin_integer_seq;
    readin_integer_seq.reserve(200);

    SeqRead input_festa(fasta_fname.c_str());
    std::ofstream reads_fhandle(reads_fname);
    uint64_t count_mapped = 0;
    int count_reads = 0;
    int inc = 0;

    for (auto festa_read: input_festa) {
        if (!(inc++ % 10))
            reads_fhandle << count_reads << endl;
        cout << festa_read->seq << endl;

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
                invalid_base_flag = true;
        }

        if (invalid_base_flag)
            continue;

        // is there a way to avoid making this copy?
        auto kmer_start_it = readin_integer_seq.begin() + readin_integer_seq.size() - kmers_size;
        auto kmer_end_it = readin_integer_seq.end();
        std::vector<uint8_t> kmer(kmer_start_it, kmer_end_it);

        if (kmer_idx.find(kmer) != kmer_idx.end() && kmer_idx_rev.find(kmer) != kmer_idx_rev.end() &&
            kmer_sites.find(kmer) != kmer_sites.end()) {
            sa_intervals = kmer_idx[kmer];
            sa_intervals_rev = kmer_idx_rev[kmer];
            sites = kmer_sites[kmer];

            it = sa_intervals.begin();
            it_rev = sa_intervals_rev.begin();

            // kmers in ref means kmers that do not cross any numbers
            // These are either in non-variable region, or are entirely within alleles
            bool first_del = !(kmers_in_ref.find(kmer) != kmers_in_ref.end());
            bool precalc_done = true;

            bidir_search_bwd(csa, (*it).first, (*it).second,
                             (*it_rev).first, (*it_rev).second,
                             readin_integer_seq.begin(),
                             readin_integer_seq.begin() + readin_integer_seq.size() - kmers_size,
                             sa_intervals, sa_intervals_rev,
                             sites, mask_allele, maxx, first_del, precalc_done);

            if (sa_intervals.size() == 1)
                //proxy for mapping is "unique horizontally"
            {
                if (!first_del)
                    // becasue sites has one element with an empty vector
                    sites.clear();

                count_mapped++;
                uint64_t no_occ = (*it).second - (*it).first;
                it = sa_intervals.begin();

                for (auto ind = (*it).first; ind < (*it).second; ind++) {
                    if (sites.empty()) {
                        if (mask_allele[csa[ind]] != 0) {
                            covgs[(mask_sites[csa[ind]] - 5) / 2][mask_allele[csa[ind]] - 1]++;
                            assert(mask_allele[csa[ind]] == mask_allele[csa[ind] + readin_integer_seq.size() - 1]);
                        }
                    } else {
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
                        if (!invalid) {
                            for (auto it_s : sites) {
                                for (auto site_pair : it_s) {
                                    auto site = site_pair.first;
                                    auto allele = site_pair.second;
                                    if (it_s != sites.back() && it_s != sites.front())
                                        assert(allele.size() == 1);
                                    // mask_a[csa[ind]] can be 0 here if the match is
                                    // coming from a skipped start_site marker
                                    if ((allele.empty()) && (mask_allele[csa[ind]] > 0))
                                        covgs[(site - 5) / 2][mask_allele[csa[ind]] - 1]++;
                                    else
                                        for (auto al : allele)
                                            covgs[(site - 5) / 2][al - 1]++;
                                }
                            }
                        }
                    }
                }
            }
            sa_intervals.clear();
            sa_intervals_rev.clear();
            sites.clear();
        }
        count_reads++;
        readin_integer_seq.clear();
    }
    reads_fhandle.close();
    cout << "Finished mapping:" << endl;
    timestamp();

    cout << count_mapped << endl;

    std::ofstream allele_coverage_fhandle(allele_coverage_fname);
    for (uint32_t i = 0; i < covgs.size(); i++) {
        for (uint32_t j = 0; j < covgs[i].size(); j++)
            allele_coverage_fhandle << covgs[i][j] << " ";
        allele_coverage_fhandle << endl;
    }
    allele_coverage_fhandle.close();

    timestamp();
    return 0;
}
