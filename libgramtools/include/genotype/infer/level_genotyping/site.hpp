#ifndef LVLGT_SITE
#define LVLGT_SITE

#include <variant>
#include "common/utils.hpp"
#include "genotype/infer/interfaces.hpp"
#include "genotype/infer/json_spec/site_spec.hpp"

namespace gram::genotype::infer {
    struct gtype_information{
        allele_vector alleles;
        GenotypeOrNull genotype;
        allele_coverages allele_covs;
        double gt_conf;
        std::size_t total_coverage;
        AlleleIds haplogroups;
    };

    class LevelGenotypedSite : public GenotypedSite {
        double gt_conf = 0.; /**< Difference in log likelihood between most likely and next most likely genotype **/
    public:
        LevelGenotypedSite() {
            auto json_site_ptr = std::make_shared<LevelGenotyped_Json_Site>();
            json_site = json_site_ptr;
        }

        ~LevelGenotypedSite() override = default;

        gtype_information get_all_gtype_info() const{
            return gtype_information{
                this->alleles,
                this->genotype,
                this->allele_covs,
                this->gt_conf,
                this->total_coverage
            };
        }

        void populate_site(gtype_information const& gtype_info){
            this->alleles = gtype_info.alleles;
            this->genotype = gtype_info.genotype;
            this->allele_covs = gtype_info.allele_covs;
            this->gt_conf = gtype_info.gt_conf;
            this->total_coverage = gtype_info.total_coverage;
            this->haplogroups = gtype_info.haplogroups;
        }

        void add_model_specific_JSON(JSON& input_json) override;
    };
}

#endif //LVLGT_SITE
