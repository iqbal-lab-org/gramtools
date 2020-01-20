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
using poisson_pmf = gram::genotype::infer::probabilities::PoissonLogPmf;
using poisson_pmf_ptr = poisson_pmf*;

namespace gram::genotype::infer {
    using numCredibleCounts = std::size_t;

    struct likelihood_related_stats{
        double mean_cov_depth, mean_pb_error, log_no_zero, log_no_zero_half_depth;
        CovCount credible_cov_t; /**< minimum coverage count to qualify as actual coverage (per-base)*/
        mutable poisson_pmf poisson_full_depth;
        mutable poisson_pmf poisson_half_depth;
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
        likelihood_related_stats const* l_stats;

        // Computed at construction time
        PerAlleleCoverage haploid_allele_coverages; /**< Coverage counts compatible with single alleles */
        PerAlleleCoverage singleton_allele_coverages; /**< Coverage counts unique to single alleles */
        std::size_t total_coverage;

        // Computed at run time
        std::multimap<double, GtypedIndices> likelihoods;
        std::shared_ptr<LevelGenotypedSite> genotyped_site; // What the class will build

    public:
        LevelGenotyperModel() : alleles(nullptr), gp_counts(nullptr) {}
        LevelGenotyperModel(allele_vector const *alleles, GroupedAlleleCounts const *gp_counts, Ploidy ploidy,
                            likelihood_related_stats const *l_stats);

        gt_site_ptr get_site() override { return std::static_pointer_cast<gt_site>(genotyped_site); }

        // Allele-level coverage
        void set_haploid_coverages(GroupedAlleleCounts const& gp_counts, AlleleId num_haplogroups);
        std::pair<float, float> compute_diploid_coverage(GroupedAlleleCounts const& gp_counts, AlleleIds ids);
        std::size_t count_total_coverage(GroupedAlleleCounts const& gp_counts);

        std::size_t count_num_haplogroups(allele_vector const& alleles);

        // Per-base coverage
        numCredibleCounts count_credible_positions(CovCount const& credible_cov_t, Allele const& allele);

        // log-likelihoods
        void compute_haploid_log_likelihoods();
        void compute_homozygous_log_likelihoods();
        void compute_heterozygous_log_likelihoods();

        // Trivial Getters
        PerAlleleCoverage const& get_haploid_covs() const {return haploid_allele_coverages;}
        PerAlleleCoverage const& get_singleton_covs() const {return singleton_allele_coverages;}
    };
}
#endif //GTYPING_MODELS
