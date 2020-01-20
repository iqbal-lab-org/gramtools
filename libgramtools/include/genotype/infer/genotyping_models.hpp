/**
 * @file Interfaces to genotyping models
 * They work on single sites: see `genotyped_site.hpp`.
 */
#ifndef GTYPING_MODELS
#define GTYPING_MODELS

#include "types.hpp"
#include "genotype/quasimap/coverage/types.hpp"
#include "genotype/infer/genotyped_site.hpp"
#include "genotype/infer/probabilities.hpp"

using namespace gram;
using poisson_pmf_ptr = gram::genotype::infer::probabilities::PoissonLogPmf*;

namespace gram::genotype::infer {
    using numCredibleCounts = std::size_t;

    struct likelihood_related_stats{
        double mean_cov_depth, mean_pb_error, log_no_zero, log_no_zero_half_depth;
        CovCount credible_cov_t; /**< minimum coverage count to qualify as actual coverage (per-base)*/
    };

    class AbstractGenotypingModel {
        virtual gt_site_ptr get_site() = 0;
    };

    /**
    Genotyping model using:
      * coverage equivalence-classes
      * alternative alleles all at same nesting level
      * genotype confidence using likelihood ratios
    */
    class LevelGenotyperModel : AbstractGenotypingModel {
        allele_vector const* alleles;
        GroupedAlleleCounts const *gp_counts;
        Ploidy ploidy;
        poisson_pmf_ptr poisson_prob;
        likelihood_related_stats const* l_stats;

        // Computed at construction time
        PerAlleleCoverage haploid_allele_coverages; /**< Coverage counts compatible with single alleles */
        PerAlleleCoverage singleton_allele_coverages; /**< Coverage counts unique to single alleles */

        // Computed at run time
        std::multimap<double, AlleleIds> likelihoods;
        std::shared_ptr<LevelGenotypedSite> genotyped_site; // What the class will build

    public:
        LevelGenotyperModel() : alleles(nullptr), gp_counts(nullptr) {}
        LevelGenotyperModel(allele_vector const* alleles, GroupedAlleleCounts const* gp_counts, Ploidy ploidy,
                            poisson_pmf_ptr poisson_prob, likelihood_related_stats const* l_stats);

        gt_site_ptr get_site() override { return std::static_pointer_cast<gt_site>(genotyped_site); }

        numCredibleCounts count_credible_positions(CovCount const& credible_cov_t, Allele const& allele);

        // Coverage
        void set_haploid_coverages(GroupedAlleleCounts const& gp_counts, AlleleId num_haplogroups);
        std::pair<float, float> compute_diploid_coverage(GroupedAlleleCounts const& gp_counts, AlleleIds ids);

        // log-likelihoods
        void haploid_log_likelihood(AlleleId const& allele);
        void diploid_log_likelihood(AlleleIds const& alleles);

        // Trivial Getters
        PerAlleleCoverage const& get_haploid_covs() const {return haploid_allele_coverages;}
        PerAlleleCoverage const& get_singleton_covs() const {return singleton_allele_coverages;}
    };
}
#endif //GTYPING_MODELS
