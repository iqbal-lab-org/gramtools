#include <genotype/infer/output_specs/segment_tracker.hpp>
#include "genotype/infer/personalised_reference.hpp"
#include "genotype/infer/interfaces.hpp"
#include "prg/coverage_graph.hpp"

namespace gram::genotype {

allele_vector get_all_alleles_to_paste(gt_site_ptr const &site, std::size_t ploidy) {
    allele_vector result(ploidy);
    auto all_site_alleles = site->get_alleles();
    GtypedIndices gts;
    if (site->is_null()) gts = GtypedIndices(ploidy, 0);
    else gts = site->get_genotype();

    if (gts.size() != ploidy) throw InconsistentPloidyException();

    for (int j{0}; j < ploidy; j++) result.at(j) = all_site_alleles.at(gts.at(j));

    return result;
}

std::size_t get_ploidy(gt_sites const& gtyped_recs){
    // If all the sites are null genotyped, will return a ploidy of one.
    std::size_t ploidy{1};
    for (auto const &site : gtyped_recs) {
        if (!site->is_null()) {
            auto gt = site->get_genotype();
            ploidy = gt.size();
            break;
        }
    }
    return ploidy;
}

// Helper function for get_personalised_ref()
void add_segment_IDs(Fastas& p_refs, std::size_t const offset, std::size_t const ploidy,
        std::string const& ID){
    for (int i{0}; i < ploidy; i++){
        std::string qualified_ID(ID + "_" + std::to_string(i));
        p_refs.at(i + offset).set_ID(qualified_ID);
    }
}

// Helper function for get_personalised_ref()
std::size_t switch_segment(Fastas &p_refs, std::size_t& offset,
        std::size_t const& ploidy, SegmentTracker &tracker) {
    if (tracker.edge() != tracker.global_edge()) {
        auto new_ID = tracker.get_ID(tracker.edge() + 1);
        offset += ploidy;
        add_segment_IDs(p_refs, offset, ploidy, new_ID);
    }
    return tracker.edge();
}

// Helper function for get_personalised_ref()
void add_invariant_sequence(Fastas& p_refs, std::size_t const offset,
                            std::size_t const ploidy, std::string const seq){
    for (int i{0}; i < ploidy; i++)
        p_refs.at(i + offset).add_sequence(seq);
}

Fastas get_personalised_ref(covG_ptr graph_root, gt_sites const &genotyped_records,
        SegmentTracker &tracker) {
    auto ploidy = get_ploidy(genotyped_records);
    auto num_segments = tracker.num_segments();
    Fastas p_refs(num_segments * ploidy);
    gram::covG_ptr cur_Node{graph_root};

    std::size_t offset{0};
    auto cur_edge = tracker.edge();
    add_segment_IDs(p_refs, offset, ploidy, tracker.get_ID(cur_edge));

    while (cur_Node->get_edges().size() > 0){
       if (cur_Node->is_bubble_start()){
           auto site_index = siteID_to_index(cur_Node->get_site_ID());
           auto site = genotyped_records.at(site_index);
           auto to_paste_alleles = get_all_alleles_to_paste(site, ploidy);
           for (int i{0}; i < ploidy; i++)
               p_refs.at(i + offset).add_sequence(to_paste_alleles.at(i).sequence);

           cur_Node = site->get_site_end_node();
           if (cur_edge == cur_Node->get_pos() - 1)
               cur_edge = switch_segment(p_refs, offset, ploidy, tracker);
       }

       if (cur_Node->has_sequence()){
           std::size_t cur_pos = cur_Node->get_pos();
           std::size_t end_pos = cur_pos + cur_Node->get_sequence_size() - 1;
           auto sequence = cur_Node->get_sequence();
           while (cur_pos <= end_pos){
                if (cur_edge <= end_pos){
                    add_invariant_sequence(p_refs, offset, ploidy,
                            sequence.substr(cur_pos - cur_Node->get_pos(), cur_edge - cur_pos + 1));
                    cur_pos = cur_edge + 1;
                    cur_edge = switch_segment(p_refs, offset, ploidy, tracker);
                }
                else {
                    add_invariant_sequence(p_refs, offset, ploidy,
                            sequence.substr(cur_pos - cur_Node->get_pos()));
                    cur_pos = end_pos + 1;
                }
            }
       }

       assert(cur_Node->get_edges().size() == 1);
       cur_Node = cur_Node->get_edges().at(0);
    }

    return p_refs;
}


    bool operator <(const Fasta& first, const Fasta& second){
    return first.sequence < second.sequence;
}

std::ostream& operator<<(std::ostream& out_stream, const Fasta& input){
   out_stream << '>' << input.ID << " " << input.desc;
    if (input.desc.back() != '\n'){
       out_stream << std::endl;
    }

    auto seq_write = input.sequence.c_str();
    auto remaining = input.sequence.size();
    while (remaining > FASTA_LWIDTH){
       out_stream.write(seq_write, FASTA_LWIDTH);
       seq_write += FASTA_LWIDTH;
       remaining -= FASTA_LWIDTH;
       out_stream << std::endl;
    }

    out_stream.write(seq_write, remaining);
    return out_stream;
}

void add_description(Fastas &p_refs, std::string const &desc) {
    int ref_number{1};
    for (auto& p_ref : p_refs) {
        p_ref.set_desc(desc);
        ref_number++;
    }
}

}
