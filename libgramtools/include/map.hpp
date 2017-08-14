#include "sequence_read/seqread.hpp"

#include "bwt_search.hpp"
#include "parameters.hpp"
#include "masks.hpp"
#include "kmers.hpp"
#include "fm_index.hpp"
#include "ranks.hpp"


#ifndef GRAMTOOLS_MAP_HPP
#define GRAMTOOLS_MAP_HPP


int map_reads(KmersData &kmers, MasksParser &masks, const Parameters &params, const FM_Index &fm_index,
              const DNA_Rank &rank_all);

std::vector<uint8_t> int_encode_read(const GenomicRead &read_sequence);

int process_read(const GenomicRead &read_sequence, int &in_sites, std::unordered_set<uint64_t> &repeats, KmersData &kmers,
                 MasksParser &masks, const Parameters &params, const DNA_Rank &rank_all, const FM_Index &fm_index);

void output_allele_coverage(Parameters &params, MasksParser &masks);


#endif //GRAMTOOLS_MAP_HPP
