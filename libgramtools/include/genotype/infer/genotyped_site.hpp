#include <prg/types.hpp>
#include "common/utils.hpp"
#include "types.hpp"

namespace gram::genotype::infer{
class AbstractGenotypedSite{
    allele_vector alleles;
    AlleleIds genotype;
    covG_ptr site_end_node;

public:
    virtual ~AbstractGenotypedSite() {};
    virtual AlleleIds const get_genotype() const = 0;
    virtual allele_vector const get_alleles() const = 0;
    virtual covG_ptr const get_site_end_node() const = 0;
};
}
