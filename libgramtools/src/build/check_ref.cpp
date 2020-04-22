#include <boost/iostreams/categories.hpp> // input_filter_tag
#include <boost/iostreams/operations.hpp> // get, WOULD_BLOCK
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp> // decompressor

#include "build/check_ref.hpp"


using namespace boost::iostreams;

struct to_upper_filter {
    typedef char              char_type;
    typedef input_filter_tag  category;

    template<typename Source>
    int get(Source& src)
    {
        int c = boost::iostreams::get(src);
        if (c == EOF) return c;
        return std::toupper((unsigned char) c);
    }
};

gram::PrgRefChecker::PrgRefChecker(std::istream &fasta_ref_handle, coverage_Graph const &cov_graph, bool gzipped) {
    filtering_istreambuf in;
    in.push(to_upper_filter());
    if (gzipped){
        in.push(gzip_decompressor());
    }
    in.push(fasta_ref_handle);

    std::string prg_first_path = get_first_prg_path(cov_graph);
    uint32_t prg_offset{0};

    std::string line, ref_seq;
    std::istream getter{&in};
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
