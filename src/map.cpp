/* TODO:
 *  *   implement boost logging
*/

#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

#include <seqread.hpp>
#include "masks.hpp"
#include "parameters.hpp"
#include "bwt_search.h"
#include "kmers.hpp"
#include "map.hpp"


uint64_t map_festa(Parameters &params, MasksParser &masks,
                   KmersData &kmers, CSA &csa) {

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
                               masks, count_reads, kmers, count_mapped, csa);

    }
    reads_fhandle.close();
    return count_mapped;
}


void process_festa_sequence(GenomicRead *festa_read, std::vector<uint8_t> &readin_integer_seq,
                            Parameters &params, MasksParser &masks, int &count_reads,
                            KmersData &kmers, uint64_t &count_mapped, CSA &csa) {

    std::cout << festa_read->seq << std::endl;
    bool invalid_base_flag = convert_festa_to_int_seq(festa_read, readin_integer_seq);
    if (invalid_base_flag)
        // TODO: should readin_integer_seq and count_reads be modified here?
        return;

    // TODO: is there a way to avoid making this copy?
    auto kmer_start_it = readin_integer_seq.begin() + readin_integer_seq.size() - params.kmers_size;
    auto kmer_end_it = readin_integer_seq.end();
    std::vector<uint8_t> kmer(kmer_start_it, kmer_end_it);

    bool kmer_within_range = kmers.index.find(kmer) != kmers.index.end()
                             && kmers.index_reverse.find(kmer) != kmers.index_reverse.end()
                             && kmers.sites.find(kmer) != kmers.sites.end();
    if (!kmer_within_range){
        count_reads++;
        readin_integer_seq.clear();
        return;
    }

    auto sa_intervals = kmers.index[kmer];
    auto sa_intervals_rev = kmers.index_reverse[kmer];
    auto sites = kmers.sites[kmer];

    auto it = sa_intervals.begin();
    auto it_rev = sa_intervals_rev.begin();

    // kmers in ref means kmers that do not cross any numbers
    // These are either in non-variable region, or are entirely within alleles
    bool first_del = !(kmers.in_reference.find(kmer) != kmers.in_reference.end());
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
                assert(masks.allele[csa[ind]]
                       == masks.allele[csa[ind] + readin_integer_seq.size() - 1]);
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
