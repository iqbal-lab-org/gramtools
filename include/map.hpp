#include "sequence_read/seqread.hpp"


#ifndef GRAMTOOLS_MAP_HPP
#define GRAMTOOLS_MAP_HPP

std::pair<int, int> map_fastaq(Parameters &params, MasksParser &masks,
                   KmersData &kmers, CSA &csa);


bool convert_fastaq_to_int_seq(GenomicRead *fastaq_read, std::vector<uint8_t> &readin_integer_seq);


void process_fastaq_sequence(GenomicRead *fastaq_read, std::vector<uint8_t> &readin_integer_seq,
        Parameters &params, MasksParser &masks, int &count_reads,
        KmersData &kmers, uint64_t &count_attempt_mapped, CSA &csa, int &in_sites, int &no_mapped,
        std::unordered_set<int> &repeats);


void output_allele_coverage(Parameters &params, MasksParser &masks);

#endif //GRAMTOOLS_MAP_HPP
