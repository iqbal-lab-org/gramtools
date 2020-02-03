#include "genotype/infer/genotyped_site.hpp"
#include "genotype/infer/allele_extracter.hpp"

using namespace gram::genotype::infer;

allele_vector gram::genotype::infer::prepend_allele(allele_vector const& original_alleles, Allele const& to_prepend){
    allele_vector result;
    result.reserve(original_alleles.size() + 1);
    result.insert(result.end(), to_prepend);
    result.insert(result.end(), original_alleles.begin(), original_alleles.end());

    return result;
}


AlleleExtracter::AlleleExtracter(covG_ptr site_start, covG_ptr site_end, gt_sites& sites)
        : genotyped_sites(&sites){
    assert(site_start->is_bubble_start());
    AlleleId haplogroup_ID{0};

    for (auto& haplogroup_start_node : site_start->get_edges()){
        allele_vector extracted_alleles = extract_alleles(haplogroup_ID, haplogroup_start_node, site_end);
        alleles.insert(alleles.end(), extracted_alleles.begin(), extracted_alleles.end());
        haplogroup_ID++;
    }
}

allele_vector AlleleExtracter::allele_combine(allele_vector const& existing, std::size_t site_index){
    // Sanity check: site_index refers to actual site
    assert( 0 <= site_index && site_index < genotyped_sites->size() );
    gt_site_ptr referent_site = genotyped_sites->at(site_index);

    allele_vector genotyped_alleles = referent_site->extract_unique_genotyped_alleles();
    allele_vector combinations(existing.size() * genotyped_alleles.size());

    std::size_t insertion_index{0};
    for (auto const& allele : existing){
        for (auto const& genotyped_allele : genotyped_alleles){
            combinations.at(insertion_index) = allele + genotyped_allele;
            ++insertion_index;
        }
    }
    return combinations;
}

void AlleleExtracter::allele_paste(allele_vector& existing, covG_ptr sequence_node){
    Allele to_paste_allele{sequence_node->get_sequence(), sequence_node->get_coverage()};
    for (auto& allele : existing) allele = allele + to_paste_allele;
}


allele_vector AlleleExtracter::extract_alleles(AlleleId const haplogroup, covG_ptr haplogroup_start, covG_ptr site_end) {
    allele_vector haplogroup_alleles{
            {"", {}, haplogroup}
    }; // Make one empty allele as starting point
    covG_ptr cur_Node{haplogroup_start};
    bool first_advance{true};

    Allele ref_allele{"", {}, 0};
    Allele to_paste_allele;

    while(cur_Node != site_end){
        if (cur_Node->is_bubble_start()) {
            auto site_index = siteID_to_index(cur_Node->get_site_ID());
            haplogroup_alleles = allele_combine(haplogroup_alleles, site_index);

            auto referent_site = genotyped_sites->at(site_index);
            cur_Node = referent_site->get_site_end_node(); // Move past site, to bubble end

            if (haplogroup == 0) { // REF extraction
                auto ref_containing_alleles = referent_site->get_alleles();
                assert(ref_containing_alleles.size() > 0);
                to_paste_allele = ref_containing_alleles.at(0);
            }
        }
        else if (! first_advance && cur_Node->is_bubble_end()) ;

        // We have an allele. Note this condition allows for a direct deletion, which will paste an empty sequence
        else {
            allele_paste(haplogroup_alleles, cur_Node);
            if (haplogroup == 0) to_paste_allele = Allele{cur_Node->get_sequence(), cur_Node->get_coverage()};
        }

        if (haplogroup == 0) ref_allele = ref_allele + to_paste_allele;

        assert(cur_Node->get_num_edges() == 1); // The only nodes with >1 neighbour are bubble starts and we skipped past those.
        cur_Node = *(cur_Node->get_edges().begin()); // Advance to the next node
        if (first_advance) first_advance = false;
    }

    if (haplogroup == 0) {
        _ref_allele_got_made_naturally = ref_allele == haplogroup_alleles.at(0);
        if (! _ref_allele_got_made_naturally )
            haplogroup_alleles = prepend_allele(haplogroup_alleles, ref_allele);
    }

    return haplogroup_alleles;
}
