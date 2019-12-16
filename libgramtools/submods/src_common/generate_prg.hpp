#ifndef COMMON_GEN_PRG
#define COMMON_GEN_PRG

#include "prg/prg_info.hpp"
#include "common/utils.hpp"
#include "kmer_index/masks.hpp"

gram::PRG_Info generate_prg_info(const marker_vec &prg_raw);

#endif //COMMON_GEN_PRG
