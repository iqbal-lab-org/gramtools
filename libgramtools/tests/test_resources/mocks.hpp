#include "common/random.hpp"
#include "gmock/gmock.h"

class MockRandomGenerator : public gram::RandomGenerator {
 public:
  MOCK_METHOD(uint32_t, generate, (uint32_t, uint32_t), (const override));
};
