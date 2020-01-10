#include "genotype/infer/genotyped_site.hpp"
#include "genotype/infer/allele_extracter.hpp"

using namespace gram::genotype::infer;

AlleleExtracter::AlleleExtracter(covG_ptr site_start, covG_ptr site_end, gt_sites& sites)
: genotyped_sites(genotyped_sites){
    ;
}

allele_vector AlleleExtracter::allele_combine(allele_vector existing, std::size_t site_index){
    // Sanity check: site_index refers to actual site
    assert( 0 <= site_index && site_index < genotyped_sites.size() );
    gtSite_ptr referent_site = genotyped_sites.at(site_index);
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
