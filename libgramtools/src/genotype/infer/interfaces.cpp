#include "genotype/infer/interfaces.hpp"

namespace gram::genotype::infer {

    /*
     * Sites
     */

allele_vector const AbstractGenotypedSite::get_unique_genotyped_alleles(allele_vector const &all_alleles,
                                                                        GenotypeOrNull const &genotype) const {

    std::set<GtypedIndex> distinct_genotypes;
    if (auto valid_gtype = std::get_if<GtypedIndices>(&genotype)) {
        // NOTE/CRUCIAL: this sorts the genotypes (eg 1,0 goes to 0, 1), which is REQUIRED for REF allele production
        distinct_genotypes = std::set<GtypedIndex>(valid_gtype->begin(), valid_gtype->end());
    } else distinct_genotypes = std::set<GtypedIndex>{0}; // If null genotype, take the reference only

    allele_vector result(distinct_genotypes.size());

    std::size_t insertion_index{0};
    for (auto const &allele_id : distinct_genotypes) {
        result.at(insertion_index) = all_alleles.at(allele_id);
        ++insertion_index;
    }
    return result;
}

AlleleIds const AbstractGenotypedSite::get_nonGenotyped_haplogroups() const{
    assert(! is_null());
    assert(alleles.size() > 0);
    assert(num_haplogroups > 0);
    AlleleIds result;

    AlleleIdSet genotyped_haplogroups;
    for (auto const& gt : std::get<GtypedIndices>(genotype)){
       genotyped_haplogroups.insert(alleles.at(gt).haplogroup);
    }

    for (AlleleId i{0}; i < num_haplogroups; i++){
        if (genotyped_haplogroups.find(i) == genotyped_haplogroups.end())
            result.push_back(i);
    }
    return result;
}
}
