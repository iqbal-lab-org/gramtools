#include <boost/iostreams/filtering_streambuf.hpp>

#include "build/check_ref.hpp"
#include "common/file_read.hpp"

using namespace boost::iostreams;

gram::PrgRefChecker::PrgRefChecker(std::istream &fasta_ref_handle, coverage_Graph const &cov_graph,
                                   bool const gzipped) {
    boost::iostreams::filtering_istreambuf in;
    input_fasta(in, fasta_ref_handle, gzipped);
    std::istream getter{&in};

    std::string prg_first_path = get_first_prg_path(cov_graph);
    uint32_t prg_offset{0};
    std::string line, ref_seq;
    while(std::getline(getter, line)){
        if (line[0] == '>') continue;
        ref_seq += line;
        auto prg_slice = prg_first_path.substr(prg_offset, line.size());
        if (prg_slice != line){
            throw std::runtime_error(
                    "Reference sequence " + line + " does not match prg slice " +
                    prg_slice + " from position " + std::to_string(prg_offset)
            );
        }
        prg_offset += line.size();
    }
    assert(prg_offset > 0);
}

std::string gram::PrgRefChecker::get_first_prg_path(coverage_Graph const &cov_graph){
    std::string path;
    auto cur_Node = cov_graph.root;
    std::size_t num_outgoing = cur_Node->get_edges().size();
    while (num_outgoing > 0){
        if (cur_Node->has_sequence()) path += cur_Node->get_sequence();
        cur_Node = cur_Node->get_edges().at(0);
        num_outgoing = cur_Node->get_edges().size();
    }
    return path;
}
