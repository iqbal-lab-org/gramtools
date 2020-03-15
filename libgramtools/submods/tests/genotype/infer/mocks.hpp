#include "genotype/infer/interfaces.hpp"
#include "genotype/infer/level_genotyping/probabilities.hpp"
#include "gmock/gmock.h"

using namespace gram;
using namespace gram::genotype::infer;

class MockGenotypedSite : public GenotypedSite{
public:
    MOCK_METHOD(site_entries,  get_model_specific_entries, (), (override));
    MOCK_METHOD(void,  null_model_specific_entries, (), (override));
};

namespace gram::genotype::infer::probabilities{
    class MockPmf : public AbstractPmf{
    public:
        MOCK_METHOD(double, compute_prob, (params const& query), (const, override));
    };
}
