#include "genotype/infer/interfaces.hpp"

namespace gram::genotype::infer {

    void GenotypedSite::populate_site(gtype_information const& gtype_info){
        this->gtype_info.alleles = gtype_info.alleles;
        this->gtype_info.genotype = gtype_info.genotype;
        this->gtype_info.allele_covs = gtype_info.allele_covs;
        this->gtype_info.total_coverage = gtype_info.total_coverage;
        this->gtype_info.haplogroups = gtype_info.haplogroups;
    }

    allele_vector const GenotypedSite::get_unique_genotyped_alleles(allele_vector const &all_alleles,
                                                                    GtypedIndices const &genotype) const {

    std::set<GtypedIndex> distinct_genotypes;
    if (! is_null()){
        // NOTE/CRUCIAL: this sorts the genotypes (eg 1,0 goes to 0, 1), which is REQUIRED for REF allele production
        distinct_genotypes.insert(genotype.begin(), genotype.end());
    } else distinct_genotypes.insert(0); // If null genotype, take the reference only

    allele_vector result(distinct_genotypes.size());

    std::size_t insertion_index{0};
    for (auto const &allele_id : distinct_genotypes) {
        result.at(insertion_index) = all_alleles.at(allele_id);
        ++insertion_index;
    }
    return result;
}


AlleleIds GenotypedSite::get_genotyped_haplogroups(allele_vector const& input_alleles,
                                                   GtypedIndices const& input_gts) const{
    AlleleIds result;
    for (auto const& gt : input_gts){
       result.push_back(input_alleles.at(gt).haplogroup);
    }
    return result;
}


}
