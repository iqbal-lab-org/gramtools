/* TODO:
 *  *   implement boost logging
*/

#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

#include "masks.hpp"
#include "parameters.hpp"
#include "bwt_search.hpp"
#include "map.hpp"
#include "bidir_search_bwd.hpp"


int map_reads(Parameters &params,
              MasksParser &masks,
              KmersData &kmers,
              const FM_Index &fm_index,
              const DNA_Rank &rank_all) {

    SeqRead reads(params.reads_fpath.c_str());
    std::ofstream reads_fhandle(params.processed_reads_fpath);
    int count_reads = 0;
    int count_mapped = 0;
    int inc = 0;
    int in_sites = 0;
    std::unordered_set<int> repeats;

    std::vector<uint8_t> readin_integer_seq;
    readin_integer_seq.reserve(200);

    for (const auto * const read: reads) {
        if (inc++ % 100000 == 0)
            reads_fhandle << count_reads << std::endl;

        process_read(*read, readin_integer_seq, params,
                     masks, count_reads, count_mapped, kmers, fm_index, in_sites,
                     repeats, rank_all);
    }
    reads_fhandle.close();
    return count_mapped;
}


void process_read(const GenomicRead &read_sequence,
                  std::vector<uint8_t> &readin_integer_seq,
                  Parameters &params,
                  MasksParser &masks,
                  int &count_reads, int &count_mapped,
                  KmersData &kmers,
                  const FM_Index &fm_index,
                  int &in_sites,
                  std::unordered_set<int> &repeats,
                  const DNA_Rank &rank_all) {

    bool invalid_base_flag = int_encode_read(read_sequence, readin_integer_seq);
    if (invalid_base_flag)
        // TODO: should readin_integer_seq and count_reads be modified here?
        return;
    auto no_occ = 0;
    std::vector<uint8_t> kmer(readin_integer_seq.begin() + readin_integer_seq.size() - params.kmers_size,
                              readin_integer_seq.end()); //is there a way to avoid making this copy?

    if (kmers.index.find(kmer) != kmers.index.end()
        && kmers.index_reverse.find(kmer) != kmers.index_reverse.end()
        && kmers.sites.find(kmer) != kmers.sites.end()) {

        auto sa_intervals = kmers.index[kmer];
        auto sa_intervals_rev = kmers.index_reverse[kmer];
        auto sites = kmers.sites[kmer];

        auto it = sa_intervals.begin();
        auto it_rev = sa_intervals_rev.begin();

        bool first_del = true;
        //kmers in ref means kmers that do not cross any numbers
        //These are either in non-variable region, or are entirely within alleles
        if (kmers.in_reference.find(kmer) != kmers.in_reference.end()) {
            //then the kmer does overlap a number, by definition.
            first_del = false;//no need to ignore first SA interval (if it was in the nonvar bit would ignore)
        }

        bool kmer_precalc_done = true;
        bidir_search_bwd(sa_intervals, sa_intervals_rev, (*it).first, (*it).second, (*it_rev).first, (*it_rev).second,
                         sites, first_del,
                         readin_integer_seq.begin(),
                         readin_integer_seq.begin() + readin_integer_seq.size() - params.kmers_size, masks.allele,
                         masks.max_alphabet_num, kmer_precalc_done,
                         rank_all, fm_index, 0);
        no_occ = 0;

        //proxy for mapping is "unique horizontally"
        if (sa_intervals.size() <= 100) {
            it = sa_intervals.begin();
            no_occ = (*it).second - (*it).first;
            count_mapped++;
            if (!first_del) {
                assert(sites.front().empty());//becasue matches are all in non variable part of PRG
                repeats.clear();
                in_sites = 0;
                for (auto ind = (*it).first; ind < (*it).second; ind++) {
                    if (masks.allele[fm_index[ind]] != 0) {
                        in_sites++;
                        if (repeats.count(masks.sites[fm_index[ind]]) == 0) repeats.insert(masks.sites[fm_index[ind]]);
                        assert(masks.allele[fm_index[ind]] ==
                               masks.allele[fm_index[ind] + readin_integer_seq.size() - 1]);
                    }
                }
            }

            auto it_ss = sites.begin();

            while (it != sa_intervals.end() && it_ss != sites.end()) {
                if (it == sa_intervals.begin() && sites.front().empty()) {
                    assert(!first_del);
                    for (auto ind = (*it).first; ind < (*it).second; ind++)
                        if (masks.allele[fm_index[ind]] != 0) {
                            masks.allele_coverage[(masks.sites[fm_index[ind]] - 5) / 2][masks.allele[fm_index[ind]] - 1] =
                                    masks.allele_coverage[(masks.sites[fm_index[ind]] - 5) / 2][masks.allele[fm_index[ind]] - 1] +
                                    1.0 / (no_occ - in_sites + repeats.size() + sa_intervals.size() - 1); //careful, might be dividing with more than we need to. size of sa_intervals is an overestimate of the number of horizontal matches, since a match that passed through 1st allele will be in a separate interval from other vertical matches from the same site
                            assert(masks.allele[fm_index[ind]] ==
                                   masks.allele[fm_index[ind] + readin_integer_seq.size() - 1]);
                        }
                } else if ((it == sa_intervals.begin() && first_del) || (it !=
                                                                                 sa_intervals.begin())) { //first_del=true - match in an interval starting with a number, all matches must be just to left of end marker
                    //if no_occ>1, two matches both starting at the end marker. If one crossed the start marker,
                    //sorina would have split into two SAs and here we are in one.
                    //so neither crosses the start marker, both start at the end. Since she only updates sites
                    //when you cross the left marker, it should be true that sites.front().back().second.size==0
                    auto it_s = *it_ss;
                    assert(!it_s.empty());

                    auto invalid = false;
                    for (auto site_pair : it_s) {
                        auto site = site_pair.first;
                        auto allele = site_pair.second;
                        if (site_pair != it_s.back() && allele.empty()) invalid = true;
                    }

                    if (!invalid) {
                        if (((*it).second - (*it).first) > 1)
                            assert(it_s.back().second.size() == 0); //vertically non-unique
                        for (auto site_pair : it_s) {
                            auto site = site_pair.first;
                            auto allele = site_pair.second;
                            if (site_pair != it_s.back() && site_pair != it_s.front()) assert(allele.size() == 1);
                            if (site_pair == it_s.back() && allele.empty()) {
                                for (auto ind = (*it).first; ind < (*it).second; ind++)
                                    if (masks.allele[fm_index[ind]] >
                                        0) { //mask_a[fm_index[ind]] can be 0 here if one match in the SA-interval is coming from a skipped start-site marker
                                        if (!first_del) { //take into account the number of matches in the reference seq
                                            assert((no_occ - in_sites + repeats.size() + sa_intervals.size() - 1) > 0);
                                            masks.allele_coverage[(site - 5) / 2][masks.allele[fm_index[ind]] - 1] =
                                                    masks.allele_coverage[(site - 5) / 2][masks.allele[fm_index[ind]] -
                                                                                          1] +
                                                    1.0 /
                                                    (no_occ - in_sites + repeats.size() + sa_intervals.size() - 1);
                                        } else
                                            masks.allele_coverage[(site - 5) / 2][masks.allele[fm_index[ind]] - 1] =
                                                    masks.allele_coverage[(site - 5) / 2][masks.allele[fm_index[ind]] -
                                                                                          1] +
                                                    1.0 / sa_intervals.size();
                                    }
                            } else if (!allele.empty()) {
                                assert(!allele.empty());
                                for (auto al:allele) {
                                    if (!first_del) {
                                        assert((no_occ - in_sites + repeats.size() + sa_intervals.size() - 1) > 0);
                                        masks.allele_coverage[(site - 5) / 2][al - 1] =
                                                masks.allele_coverage[(site - 5) / 2][al - 1] +
                                                1.0 / (no_occ - in_sites + repeats.size() + sa_intervals.size() - 1);
                                    } else
                                        masks.allele_coverage[(site - 5) / 2][al - 1] =
                                                masks.allele_coverage[(site - 5) / 2][al - 1] +
                                                1.0 / sa_intervals.size();
                                }
                            }
                        }
                    }
                }
                ++it;
                ++it_ss;
            }
        } else {
            no_occ = 0;
        }
        sa_intervals.clear();
        sa_intervals_rev.clear();
        sites.clear();
    } else {
        no_occ = 0;
    }

    count_reads++;
    readin_integer_seq.clear();
}


bool int_encode_read(const GenomicRead &read_sequence, std::vector<uint8_t> &readin_integer_seq) {
    bool invalid_base_flag = false;
    const auto sequence_length = strlen(read_sequence.seq);
    for (int i = 0; i < sequence_length; i++) {
        if (read_sequence.seq[i] == 'A' or read_sequence.seq[i] == 'a')
            readin_integer_seq.push_back(1);
        else if (read_sequence.seq[i] == 'C' or read_sequence.seq[i] == 'c')
            readin_integer_seq.push_back(2);
        else if (read_sequence.seq[i] == 'G' or read_sequence.seq[i] == 'g')
            readin_integer_seq.push_back(3);
        else if (read_sequence.seq[i] == 'T' or read_sequence.seq[i] == 't')
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
