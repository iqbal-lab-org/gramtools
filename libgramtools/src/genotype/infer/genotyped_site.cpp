#include "genotype/infer/genotyped_site.hpp"

namespace gram::genotype::infer {

allele_vector const get_unique_genotyped_alleles(allele_vector const &all_alleles, GenotypeOrNull const &genotype) {

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
}
