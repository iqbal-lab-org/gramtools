#include "simulate/induce_genotypes.hpp"

using namespace gram::simulate;

void NodeThread::visit(nt_ptr_v & to_visit, std::string const& sequence) const{
    auto node_size = prg_node->get_sequence_size();
    if (prg_node->has_sequence()){
        auto seq_slice = sequence.substr(offset, node_size);
        // Only add new nodes to visit if this node's sequence matches
        if (seq_slice != prg_node->get_sequence()) return;
    }
    for (auto const& n : prg_node->get_edges()){
        auto new_node = std::make_shared<const NodeThread>(shared_from_this(), n, offset + node_size);
        to_visit.push_back(new_node);
    }
}

nt_ptr gram::simulate::thread_sequence(covG_ptr root, std::string const& sequence){
    auto cur_Node = std::make_shared<const NodeThread>(nullptr, root, 0);
    nt_ptr_v to_visit{cur_Node}, endpoints;
    while (not to_visit.empty()){
        cur_Node = to_visit.back();
        to_visit.pop_back();
        if (not cur_Node->has_next())
            endpoints.push_back(cur_Node);
        else cur_Node->visit(to_visit, sequence);
    }

    std::string msg{};
    switch(endpoints.size()){
        case 0:
            msg = "Could not thread a path through the prg for sequence: \n" + sequence;
            throw NoEndpoints(msg);
        case 1:
            return endpoints.back();
        default:
            msg = "Found more than one path through the prg for sequence: \n" + sequence;
            throw TooManyEndpoints(msg);
    }
}
