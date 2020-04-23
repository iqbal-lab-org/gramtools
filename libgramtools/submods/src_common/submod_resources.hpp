#ifndef COMMON_GEN_PRG
#define COMMON_GEN_PRG

#include "prg/prg_info.hpp"
#include "common/data_types.hpp"


namespace gram::submods{
    gram::PRG_Info generate_prg_info(const marker_vec &prg_raw);

    std::string decode(const uint64_t base);
}

using namespace gram::submods;

#endif //COMMON_GEN_PRG
