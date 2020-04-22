/**
 * Takes in a fasta reference and a coverage_Graph, and checks the first path
 * in the graph corresponds to the reference.
 */
#ifndef GRAMTOOLS_CHECK_REF_HPP
#define GRAMTOOLS_CHECK_REF_HPP

#include <iostream>

#include "prg/coverage_graph.hpp"

namespace gram {
    class PrgRefChecker {
    public:
        PrgRefChecker(std::istream &fasta_ref_handle, coverage_Graph const &cov_graph,
                bool gzipped = false);

        static std::string get_first_prg_path(coverage_Graph const &cov_graph);
    };
}


#endif //GRAMTOOLS_CHECK_REF_HPP
