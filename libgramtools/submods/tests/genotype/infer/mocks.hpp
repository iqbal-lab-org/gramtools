#include "genotype/infer/genotyped_site.hpp"
#include "gmock/gmock.h"

using namespace gram;
using namespace gram::genotype::infer;

class MockGenotypedSite : public AbstractGenotypedSite{
public:
    MOCK_METHOD(AlleleIds const, get_genotype, (), (const, override));
    MOCK_METHOD(allele_vector const, get_alleles, (), (const, override));
};