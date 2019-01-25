/** @file
 * Sets up the data structures supporting vBWT search, and builds the kmer index.
 * For data structure specification see `gram::PRG_Info`.
 */
#ifndef GRAMTOOLS_BUILD_HPP
#define GRAMTOOLS_BUILD_HPP

namespace gram::commands::build {
    void run(const Parameters &parameters);
}

#endif //GRAMTOOLS_BUILD_HPP
