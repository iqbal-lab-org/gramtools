#ifndef GRAMTOOLS_SIMULATE_HPP
#define GRAMTOOLS_SIMULATE_HPP

#include "parameters.hpp"
#include "genotype/infer/interfaces.hpp"
#include "common/random.hpp"
#include "genotype/infer/output_specs/fields.hpp"

using namespace gram::genotype::infer;

namespace gram::simulate {
    class RandomGenotypedSite : public GenotypedSite {
    public:
        RandomGenotypedSite() = default;
        site_entries get_model_specific_entries() override{ return {}; }
        void null_model_specific_entries() override {}
    };

    using Seed = uint32_t;

    class RandomGenotyper : public Genotyper {
    public:
        /**
         * The genotyping process is the same in form to `gram::genotype::infer::LevelGenotyper`
         * except that genotype is randomly assigned among the list of alleles.
         */
        RandomGenotyper(coverage_Graph const &cov_graph);

        header_vec get_model_specific_headers() override;
    };

    /**
     * @param use_ref_allele : if the ref allele is not consistent with the inner bubbles,
        it is still produced by the `AlleleExtracter`. In that case this param should be false.
     */
    gt_site_ptr make_randomly_genotyped_site(RandomGenerator const *const rand, allele_vector const &alleles,
                                             bool const use_ref_allele = true);
}

namespace gram::commands::simulate {
    void run(SimulateParams const &parameters);
}

#endif //GRAMTOOLS_SIMULATE_HPP
