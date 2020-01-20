/**
 * @file Interfaces to genotyped site classes
 */
#ifndef GENOTYPED_SITE
#define GENOTYPED_SITE

#include "common/utils.hpp"
#include "types.hpp"

namespace gram::genotype::infer{
class AbstractGenotypedSite{
protected:
    allele_vector alleles;
    AlleleIds genotype;
    covG_ptr site_end_node;

public:
    virtual ~AbstractGenotypedSite() {};
    virtual AlleleIds const get_genotype() const = 0;
    virtual allele_vector const get_alleles() const = 0;
    virtual covG_ptr const get_site_end_node() const = 0;
};

class LevelGenotypedSite : public AbstractGenotypedSite{
public:
    //~LevelGenotypedSite() = default;
    AlleleIds const get_genotype() const override {return genotype;}
    allele_vector const get_alleles() const override {return alleles;}
    covG_ptr const get_site_end_node() const override {return site_end_node;}
};
}

#endif //GENOTYPED_SITE
