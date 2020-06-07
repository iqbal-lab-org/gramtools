#include "genotype/infer/allele_extracter.hpp"

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


gt_sites gram::simulate::make_nulled_sites(coverage_Graph const& input_prg) {
    using namespace gram::genotype::infer;
    gt_sites genotyped_records(input_prg.bubble_map.size());

    // Genotype each bubble in the PRG, in most nested to less nested order.
    for (auto const &bubble_pair : input_prg.bubble_map) {
        auto site = std::make_shared<SimulatedSite>();
        auto extracter = AlleleExtracter(bubble_pair.first, bubble_pair.second, genotyped_records);
        site->set_alleles({extracter.get_alleles().at(0)});
        site->set_pos(bubble_pair.first->get_pos());
        site->make_null();
        site->set_site_end_node(bubble_pair.second);

        auto site_ID = bubble_pair.first->get_site_ID();
        auto site_index = siteID_to_index(site_ID);
        genotyped_records.at(site_index) = site;
    }
    return genotyped_records;
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

Allele extract_allele(nt_ptr const& end_point, Marker const& target_site_ID){
    auto cur_Node = end_point->get_parent();
    auto cur_prg_node = cur_Node->get_prg_node();
    std::string sequence;
    AlleleId haplogroup{ALLELE_UNKNOWN};
    while (true){
        if (cur_prg_node->is_bubble_start() && cur_prg_node->get_site_ID() == target_site_ID)
            break;
        if (haplogroup == ALLELE_UNKNOWN && cur_prg_node->get_site_ID() == target_site_ID)
            haplogroup = cur_prg_node->get_allele_ID();
        sequence.insert(0, cur_prg_node->get_sequence());
        cur_Node = cur_Node->get_parent();
        cur_prg_node = cur_Node->get_prg_node();
    }
    return Allele{sequence, {}, haplogroup};
}


void gram::simulate::apply_genotypes(nt_ptr const& end_point, gt_sites const& sites){
    auto cur_Node = end_point;
    auto cur_prg_node = cur_Node->get_prg_node();
    while (cur_Node->get_parent() != nullptr){
        if (cur_prg_node->is_bubble_end()){
            auto site_ID = cur_prg_node->get_site_ID();
            auto& site = sites.at(siteID_to_index(site_ID));
            auto extracted_allele = extract_allele(cur_Node, site_ID);
            allele_vector site_alleles{site->get_alleles()};
            if (extracted_allele.sequence == site_alleles.at(0).sequence)
              site->populate_site(gtype_information{
                 site_alleles,
                 GtypedIndices{0},
                 {0},
                 1,
                 AlleleIds{0}});
            else {
                site_alleles.push_back(extracted_allele);
                site->populate_site(gtype_information{
                        site_alleles,
                        GtypedIndices{1},
                        {0, 1},
                        1,
                        AlleleIds{0, extracted_allele.haplogroup}
                });
            }
        }
        cur_Node = cur_Node->get_parent();
        cur_prg_node = cur_Node->get_prg_node();
    }
}
