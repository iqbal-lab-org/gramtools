#ifndef LVLGT_SITE
#define LVLGT_SITE

#include <variant>
#include "common/utils.hpp"
#include "genotype/infer/interfaces.hpp"

namespace gram::genotype::infer {
    class LevelGenotypedSite : public AbstractGenotypedSite {
        double gt_conf; /**< Difference in log likelihood between most likely and next most likely genotype **/
    public:
        LevelGenotypedSite() {}

        ~LevelGenotypedSite() override = default;

        GenotypeOrNull const get_genotype() const override { return genotype; }

        allele_vector const get_alleles() const override { return alleles; }

        covG_ptr const get_site_end_node() const override { return site_end_node; }

        void make_null() override {
            genotype = false;
            gt_conf = 0.;
        }


        void set_genotype(GtypedIndices const indices, double const gt_confidence) {
            genotype = indices;
            gt_conf = gt_confidence;
        };

        void set_alleles(allele_vector const &chosen_alleles) {
            alleles = chosen_alleles;
        }


        /** Whether the site is null genotyped */
        bool is_null() const override {
            if (std::holds_alternative<bool>(genotype)) return true;
            else return false;
        };
    };
}

#endif //LVLGT_SITE
