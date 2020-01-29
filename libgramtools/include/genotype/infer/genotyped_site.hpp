/**
 * @file Interfaces to genotyped site classes
 */
#ifndef GENOTYPED_SITE
#define GENOTYPED_SITE

#include <variant>
#include "common/utils.hpp"
#include "types.hpp"

namespace gram::genotype::infer{
    using GtypedIndex = std::size_t; /**< The index of an allele in an allele vector */
    using GtypedIndices = std::vector<GtypedIndex>;
    using GenotypeOrNull = std::variant<GtypedIndices, bool>;


class AbstractGenotypedSite{
protected:
    allele_vector alleles;
    GenotypeOrNull genotype;
    covG_ptr site_end_node;
    std::size_t num_haplogroups = 0; /**< The number of outgoing edges from the bubble start */

public:
    virtual ~AbstractGenotypedSite() {};
    virtual GenotypeOrNull const get_genotype() const = 0;
    virtual allele_vector const get_alleles() const = 0;
    virtual covG_ptr const get_site_end_node() const = 0;
    virtual bool is_null() const = 0;

    void set_site_end_node(covG_ptr const& end_node) {site_end_node = end_node;}
    std::size_t const& get_num_haplogroups() {return num_haplogroups;}
    void set_num_haplogroups(std::size_t const& num_haps) {num_haplogroups = num_haps;}
    bool const has_alleles() const { return alleles.size() > 0 ;}


    /**
     * Given alleles and GT, return the alleles referred to by GT
     */
    allele_vector const get_unique_genotyped_alleles
            (allele_vector const& all_alleles, GenotypeOrNull const& genotype) const;

    /**
     * Return site's called alleles directly, from its own alleles
     */
    allele_vector const get_unique_genotyped_alleles() const{
        return get_unique_genotyped_alleles(alleles, genotype);
    }

    /**
     * Return site's called alleles indirectly
     * This version exists for allowing to use mocked alleles and genotypes
     */
    allele_vector const extract_unique_genotyped_alleles() const{
        auto extracted_alleles = this->get_alleles();
        auto extracted_gts = this->get_genotype();
        return get_unique_genotyped_alleles(extracted_alleles, extracted_gts);
    }

    /**
     * Produce the haplogroups that have not been genotyped, for use in nested
     * site invalidation.
     */
    AlleleIds const get_nonGenotyped_haplogroups() const;
};

class LevelGenotypedSite : public AbstractGenotypedSite{
    double gt_conf; /**< Difference in log likelihood between most likely and next most likely genotype **/
public:
    LevelGenotypedSite() {}
    ~LevelGenotypedSite() = default;
    GenotypeOrNull const get_genotype() const override {return genotype;}
    allele_vector const get_alleles() const override {return alleles;}
    covG_ptr const get_site_end_node() const override {return site_end_node;}


    void set_genotype(GtypedIndices const indices, double const gt_confidence){
        genotype = indices;
        gt_conf = gt_confidence;
    };

    void set_alleles(allele_vector const& chosen_alleles){
        alleles = chosen_alleles;
    }


    void make_null() {
        genotype = false;
        gt_conf = 0.;
    }

    /** Whether the site is null genotyped */
    bool is_null() const override {
        if ( std::holds_alternative<bool>(genotype) ) return true;
        else return false;
    };
};
}

#endif //GENOTYPED_SITE
