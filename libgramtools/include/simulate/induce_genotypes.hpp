#ifndef GMTOOLS_SIMU_INDUCE_GTS_HPP
#define GMTOOLS_SIMU_INDUCE_GTS_HPP

#include "genotype/infer/types.hpp"
#include "prg/coverage_graph.hpp"

#include <memory>
#include <vector>

namespace gram::simulate{
    class NodeThread;
    using nt_ptr = std::shared_ptr<const NodeThread>;
    using nt_ptr_v = std::vector<nt_ptr>;
    using gt_sites = ::genotype::infer::gt_sites;

    class NoEndpoints : public std::runtime_error{
        using std::runtime_error::runtime_error;
    };
    using TooManyEndpoints = NoEndpoints;


class NodeThread : public std::enable_shared_from_this<NodeThread const>{
    public:
        explicit NodeThread(nt_ptr const& input_parent, covG_ptr input_prg_node, int input_offset) :
            parent(input_parent), prg_node(input_prg_node), offset(input_offset) {}

        ~NodeThread() = default;

        NodeThread(NodeThread const& other) = delete;
        NodeThread(NodeThread const&& other) = delete;
        NodeThread& operator=(NodeThread const& other) = delete;
        NodeThread& operator=(NodeThread const&& other) = delete;

        covG_ptr const& get_prg_node() const{return prg_node;}
        int get_offset() const{return offset;}
        bool has_next() const {return prg_node->get_num_edges() > 0;}
        void visit(nt_ptr_v & to_visit, std::string const& sequence) const;

    private:
        nt_ptr parent;
        covG_ptr prg_node;
        int offset;
    };

    nt_ptr thread_sequence(covG_ptr root, std::string const& sequence);
    void apply_genotypes(nt_ptr end_point, gt_sites const& sites);

}

#endif //GMTOOLS_SIMU_INDUCE_GTS_HPP
