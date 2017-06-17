#include "sequence_read/seqread.hpp"

#include "bwt_search.h"
#include "parameters.hpp"
#include "masks.hpp"
#include "kmers.hpp"
#include "process_prg.hpp"


#ifndef GRAMTOOLS_MAP_HPP
#define GRAMTOOLS_MAP_HPP


uint64_t map_festa(Parameters &params, MasksParser &masks,
                   KmersData &kmers, CSA &fm_index, const VariantMarkers &variants,
                   std::unordered_map<uint8_t,std::vector<uint64_t>>& rank_all);

bool convert_festa_to_int_seq(GenomicRead *festa_read, std::vector<uint8_t> &readin_integer_seq);

void process_festa_sequence(GenomicRead *festa_read, std::vector<uint8_t> &readin_integer_seq,
                            Parameters &params, MasksParser &masks, int &count_reads,
                            KmersData &kmers, uint64_t &count_mapped, CSA &csa, int &in_sites, int &no_mapped,
                            std::unordered_set<int> &repeats, const VariantMarkers &variants,
                            std::unordered_map<uint8_t,std::vector<uint64_t>>& rank_all);

void output_allele_coverage(Parameters &params, MasksParser &masks);

#endif //GRAMTOOLS_MAP_HPP
