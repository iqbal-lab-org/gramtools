#ifndef LVLGT_MODEL
#define LVLGT_MODEL

#include "genotype/quasimap/coverage/types.hpp"
#include "site.hpp"
#include "probabilities.hpp"

using namespace gram;
using poisson_pmf = gram::genotype::infer::probabilities::PoissonLogPmf;
using poisson_pmf_ptr = poisson_pmf*;

namespace gram::genotype::infer {
    using numCredibleCounts = std::size_t;
    using multiplicities = std::vector<bool>;
    using likelihood_map = std::multimap<double, GtypedIndices, std::greater<double>>;
    using memoised_coverages = std::map<AlleleIds, allele_coverages>;

    struct likelihood_related_stats {
        double mean_cov_depth,
                mean_pb_error, log_mean_pb_error,
                log_no_zero, log_no_zero_half_depth;
        CovCount credible_cov_t; /**< minimum coverage count to qualify as actual coverage (per-base)*/
        mutable poisson_pmf poisson_full_depth;
        mutable poisson_pmf poisson_half_depth;
    };

    /**
    Genotyping model using:
      * coverage equivalence-classes
      * alternative alleles all at same nesting level
      * genotype confidence using likelihood ratios
    */
    class LevelGenotyperModel : AbstractGenotypingModel {
        allele_vector *alleles;
        GroupedAlleleCounts const *gp_counts;
        Ploidy ploidy;
        likelihood_related_stats const *l_stats;

        // Computed at construction time
        PerAlleleCoverage haploid_allele_coverages; /**< Coverage counts compatible with single alleles */
        PerAlleleCoverage singleton_allele_coverages; /**< Coverage counts unique to single alleles */
        memoised_coverages computed_coverages;
        std::size_t total_coverage;

        // Computed at run time
        likelihood_map likelihoods; // Store highest likelihoods first
        std::shared_ptr <LevelGenotypedSite> genotyped_site; // What the class will build

    public:
        LevelGenotyperModel() : gp_counts(nullptr) {}

        /**
         *
         * @param ignore_ref_allele if true, the ref allele was not produced naturally,
         * and we do not consider it for genotyping
         */
        LevelGenotyperModel(allele_vector const *input_alleles, GroupedAlleleCounts const *gp_counts,
                            Ploidy ploidy, likelihood_related_stats const *l_stats,
                            bool ignore_ref_allele = false);

        gt_site_ptr get_site() override { return std::static_pointer_cast<gt_site>(genotyped_site); }

        // Allele-level coverage
        void set_haploid_coverages(GroupedAlleleCounts const &gp_counts, AlleleId num_haplogroups);

        /**
         * Alleles with no sequence correspond to direct deletions.
         * In this case they get assigned coverage by this function,
         * using the grouped allele coverages, as if they had a single base.
         */
        void assign_coverage_to_empty_alleles(allele_vector &alleles);

        /**
         *
         * Note: Due to nesting, the alleles can be from the same haplogroup; in which case, they have the same
         * haploid coverage, and they get assigned half of it each.
         */
        std::pair<double, double> compute_diploid_coverage(GroupedAlleleCounts const &gp_counts, AlleleIds ids,
                                                           multiplicities const &haplogroup_multiplicities);

        std::size_t count_total_coverage(GroupedAlleleCounts const &gp_counts);

        // Counting
        AlleleIds get_haplogroups(allele_vector const &alleles, GtypedIndices const &gtype) const;
        std::vector<bool> count_num_haplogroups(allele_vector const &alleles);

        /**
         * Express genotypes as relative to chosen alleles.
         * For eg, {0, 2, 4} in the original set of possible alleles goes to {0, 1, 2} in the 3 called alleles (yes,
         * this is Triploid example).
         */
        GtypedIndices rescale_genotypes(GtypedIndices const &genotypes);

        // Permutations
        /**
         * Credit: https://stackoverflow.com/a/9430993/12519542
         */
        std::vector <GtypedIndices> get_permutations(const GtypedIndices &indices, std::size_t const subset_size);

        // Per-base coverage
        numCredibleCounts count_credible_positions(CovCount const &credible_cov_t, Allele const &allele);

        /**
         * Haploid genotype likelihood
         */
        void compute_haploid_log_likelihoods();

        /**
         * Diploid homozygous
         */
        void compute_homozygous_log_likelihoods(multiplicities const &haplogroup_multiplicities);

        /**
         * Diploid. Because of the large possible number of diploid combinations,
         * (eg for 10 alleles, 45), we only consider for combination those alleles
         * that have at least one unit of coverage unique to them.
         */
        void compute_heterozygous_log_likelihoods(multiplicities const &haplogroup_multiplicities);

        // Trivial Getters
        PerAlleleCoverage const &get_haploid_covs() const { return haploid_allele_coverages; }

        PerAlleleCoverage const &get_singleton_covs() const { return singleton_allele_coverages; }

        likelihood_map const &get_likelihoods() const { return likelihoods; }
    };
}

#endif //LVLGT_MODEL
