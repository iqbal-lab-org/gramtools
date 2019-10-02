/** @file
 * Sets up the data structures supporting vBWT search, and builds the kmer index.
 * For data structure specification see `gram::PRG_Info`.
 */

#ifndef GRAMTOOLS_BUILD_HPP
#define GRAMTOOLS_BUILD_HPP

#include "common/parameters.hpp"
#include "common/timer_report.hpp"

#include "prg/prg.hpp"
#include "prg/load_PRG_string.hpp"
#include "prg/make_data_structures.hpp"
#include "prg/masks.hpp"

#include "kmer_index/build.hpp"
#include "kmer_index/dump.hpp"

namespace gram::commands::build {
    void run(const Parameters &parameters);
}

#endif //GRAMTOOLS_BUILD_HPP
