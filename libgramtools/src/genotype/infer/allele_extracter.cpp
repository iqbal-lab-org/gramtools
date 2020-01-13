#include "genotype/infer/genotyped_site.hpp"
#include "genotype/infer/allele_extracter.hpp"

using namespace gram::genotype::infer;

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
    gtSite_ptr referent_site = genotyped_sites->at(site_index);
    auto referent_genotype = referent_site->get_genotype();
    auto referent_alleles = referent_site->get_alleles();

    std::set<AlleleId> distinct_genotypes(referent_genotype.begin(), referent_genotype.end());
    allele_vector combinations(existing.size() * distinct_genotypes.size());

    std::size_t insertion_index{0};
    for (auto const& allele : existing){
        for (auto const& allele_id : distinct_genotypes){
            combinations.at(insertion_index) = allele + referent_alleles.at(allele_id);
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

    while(cur_Node != site_end){
        if (cur_Node->is_bubble_start()) {
            auto site_index = siteID_to_index(cur_Node->get_site_ID());
            haplogroup_alleles = allele_combine(haplogroup_alleles, site_index);

            auto referent_site = genotyped_sites->at(site_index);
            cur_Node = referent_site->get_site_end_node(); // Move past site, to bubble end
        }
        else if (cur_Node->has_sequence()) allele_paste(haplogroup_alleles, cur_Node);

        assert(cur_Node->get_num_edges() == 1); // The only nodes with >1 neighbour are bubble starts and we skipped past those.
        cur_Node = *(cur_Node->get_edges().begin()); // Advance to the next node
    }

    return haplogroup_alleles;
}
