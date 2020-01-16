/**
 * @file Interfaces to genotyping models
 * They work on single sites: see `genotyped_site.hpp`.
 */
#ifndef LVL_GTYPER
#define LVL_GTYPER

#include "types.hpp"
#include "genotype/quasimap/coverage/types.hpp"
#include "genotype/infer/genotyped_site.hpp"

using namespace gram;

namespace gram::genotype::infer {
    using numCredibleCounts = std::size_t;

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
        GroupedAlleleCounts const *gp_counts;
        Ploidy ploidy;
        CovCount credible_cov_t; /**< minimum coverage count to qualify as actual coverage*/

        PerAlleleCoverage haploid_allele_coverages; /**< Coverage counts compatible with single alleles */
        PerAlleleCoverage singleton_allele_coverages; /**< Coverage counts unique to single alleles */

    public:
        LevelGenotyper() : gp_counts(nullptr) {}
        gt_site_ptr get_site() override { return std::static_pointer_cast<gt_site>(genotyped_site); }

        void set_haploid_coverages(GroupedAlleleCounts const& gp_counts, AlleleId num_haplogroups);
        std::pair<float, float> compute_diploid_coverage(GroupedAlleleCounts const& gp_counts, AlleleIds ids);
        numCredibleCounts count_credible_positions(CovCount const& credible_cov_t, Allele const& allele);

        // Trivial Getters
        PerAlleleCoverage const& get_haploid_covs() const {return haploid_allele_coverages;}
        PerAlleleCoverage const& get_singleton_covs() const {return singleton_allele_coverages;}
    };
}

#endif //LVL_GTYPER
