#include "sequence_read/seqread.hpp"


#ifndef GRAMTOOLS_MAP_HPP
#define GRAMTOOLS_MAP_HPP

uint64_t map_festa(Parameters &params, MasksParser &masks,
                   KmersData &kmers, CSA &csa);


bool convert_festa_to_int_seq(GenomicRead *festa_read, std::vector<uint8_t> &readin_integer_seq);


void process_festa_sequence(GenomicRead *festa_read, std::vector<uint8_t> &readin_integer_seq,
        Parameters &params, MasksParser &masks, int &count_reads,
        KmersData &kmers, uint64_t &count_attempt_mapped, CSA &csa, int &in_sites, int &no_mapped,
        std::unordered_set<int> &repeats);


void output_allele_coverage(Parameters &params, MasksParser &masks);

#endif //GRAMTOOLS_MAP_HPP
