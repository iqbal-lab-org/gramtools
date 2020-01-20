/**
 * @file Interfaces to genotyped site classes
 */
#ifndef GENOTYPED_SITE
#define GENOTYPED_SITE

#include <variant>
#include "common/utils.hpp"
#include "types.hpp"

namespace gram::genotype::infer{
    using GtypedIndex = std::size_t;
using GtypedIndices = std::vector<GtypedIndex>;
using GenotypeOrNull = std::variant<GtypedIndices, bool>;

class AbstractGenotypedSite{
protected:
    allele_vector alleles;
    GenotypeOrNull genotype;
    covG_ptr site_end_node;

public:
    virtual ~AbstractGenotypedSite() {};
    virtual GenotypeOrNull const get_genotype() const = 0;
    virtual allele_vector const get_alleles() const = 0;
    virtual covG_ptr const get_site_end_node() const = 0;

};

class LevelGenotypedSite : public AbstractGenotypedSite{
    double gt_conf; /**< Difference in log likelihood between most likely and next most likely genotype **/
public:
    //~LevelGenotypedSite() = default;
    GenotypeOrNull const get_genotype() const override {return genotype;}
    allele_vector const get_alleles() const override {return alleles;}
    covG_ptr const get_site_end_node() const override {return site_end_node;}

    void set_genotype(GtypedIndices const indices, double const gt_confidence){
        genotype = indices;
        gt_conf = gt_confidence;
    };


    void make_null() {
        genotype = false;
        gt_conf = 0.;
    }

    /** Whether the site is null genotyped */
    bool is_null() const {
        if ( std::holds_alternative<bool>(genotype) ) return true;
        else return false;
    };
};
}

#endif //GENOTYPED_SITE
