#include "sequence_read/seqread.hpp"

#include "bwt_search.hpp"
#include "parameters.hpp"
#include "masks.hpp"
#include "kmers.hpp"
#include "fm_index.hpp"
#include "ranks.hpp"


#ifndef GRAMTOOLS_MAP_HPP
#define GRAMTOOLS_MAP_HPP


uint64_t map_reads(Parameters &params, MasksParser &masks,
                   KmersData &kmers, const FM_Index &fm_index,
                   const DNA_Rank &rank_all);

bool int_encode_read(GenomicRead *read_sequence, std::vector<uint8_t> &readin_integer_seq);

void process_read(GenomicRead *read_sequence, std::vector<uint8_t> &readin_integer_seq,
                  Parameters &params, MasksParser &masks, int &count_reads, int &count_mapped,
                  KmersData &kmers, const FM_Index &fm_index, int &in_sites, std::unordered_set<int> &repeats,
                  const DNA_Rank &rank_all);

void output_allele_coverage(Parameters &params, MasksParser &masks);

#endif //GRAMTOOLS_MAP_HPP
