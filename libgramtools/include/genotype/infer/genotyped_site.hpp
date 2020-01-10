#include "common/utils.hpp"
#include "types.hpp"

class AbstractGenotypedSite{
protected:
    virtual ~AbstractGenotypedSite();
    virtual allele_vector getGenotypedAlleles();
    allele_vector alleles;
    std::vector<AlleleId> genotype;

public:
    std::vector<AlleleId> const& get_genotype() const{
        return genotype;
    }
    allele_vector const& get_alleles() const{
        return alleles;
    }
};
}
