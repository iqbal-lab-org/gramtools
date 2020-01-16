/**
 * @file Interfaces to genotyping models
 */
#ifndef LVL_GTYPER
#define LVL_GTYPER

#include "types.hpp"
#include "genotype/quasimap/coverage/types.hpp"
#include "genotype/infer/genotyped_site.hpp"

using namespace gram;

namespace gram::genotype::infer {

    class AbstractGenotypingModel {
        virtual gt_site_ptr get_site() = 0;
    };

    /**
    Genotyping model using:
      * coverage equivalence-classes
      * alternative alleles all at same nesting level
      * genotype confidence using likelihood ratios
    */
    class LevelGenotyper : AbstractGenotypingModel {
        std::shared_ptr<LevelGenotypedSite> genotyped_site; // Build with make_shared
        allele_vector alleles;
        Coverage const *coverage_struct;
        Ploidy ploidy;

        PerAlleleCoverage haploid_allele_coverages; /**< Coverage counts compatible with single alleles */
        PerAlleleCoverage singleton_allele_coverages; /**< Coverage counts unique to single alleles */

    public:
        LevelGenotyper() : coverage_struct(nullptr) {}
        gt_site_ptr get_site() override { return std::static_pointer_cast<gt_site>(genotyped_site); }

        void set_haploid_coverages(Coverage const * coverage_struct);
    };
}

#endif //LVL_GTYPER
