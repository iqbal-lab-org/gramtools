#ifndef LVLGT_SITE
#define LVLGT_SITE

#include <variant>
#include "common/utils.hpp"
#include "genotype/infer/interfaces.hpp"

namespace gram::genotype::infer {
    struct gtype_information{
        allele_vector alleles;
        GenotypeOrNull genotype;
        allele_coverages allele_covs;
        double gt_conf;
        std::size_t total_coverage;
    };

    class LevelGenotypedSite : public AbstractGenotypedSite {
        allele_coverages allele_covs;
        double gt_conf; /**< Difference in log likelihood between most likely and next most likely genotype **/
        std::size_t total_coverage; /**< Total coverage on this site */
    public:
        LevelGenotypedSite() {}

        ~LevelGenotypedSite() override = default;

        GenotypeOrNull const get_genotype() const override { return genotype; }
        allele_vector const get_alleles() const override { return alleles; }
        covG_ptr const get_site_end_node() const override { return site_end_node; }
        gtype_information get_all_gtype_info() const{
            return gtype_information{
                this->alleles,
                this->genotype,
                this->allele_covs,
                this->gt_conf,
                this->total_coverage
            };
        }

        void make_null() override {
            this->genotype = false;
            this->gt_conf = 0.;
        }

        void populate_site(gtype_information const& gtype_info){
            this->alleles = gtype_info.alleles;
            this->genotype = gtype_info.genotype;
            this->allele_covs = gtype_info.allele_covs;
            this->gt_conf = gtype_info.gt_conf;
            this->total_coverage = gtype_info.total_coverage;
        }

        void set_alleles(allele_vector const& alleles){ this->alleles = alleles; };
        void set_genotype(GtypedIndices const& gtype, double gt_conf){
            this->genotype = gtype;
            this->gt_conf = gt_conf;
        };
        void set_total_coverage(std::size_t const& total_cov){total_coverage = total_cov;}

        /** Whether the site is null genotyped */
        bool is_null() const override {
            if (std::holds_alternative<bool>(genotype)) return true;
            else return false;
        };
    };
}

#endif //LVLGT_SITE
