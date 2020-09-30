#ifndef GRAMTOOLS_RANDOM_HPP
#define GRAMTOOLS_RANDOM_HPP

#include <cstdint>
#include <random>

#include "genotype/parameters.hpp"

namespace gram {

/**
 * Abstract base class used for mocking in unit tests
 */
class RandomGenerator {
 public:
  virtual ~RandomGenerator(){};

  virtual uint32_t generate(uint32_t min, uint32_t max) = 0;

  SeedSize operator()() { return random_number_generator(); }

 protected:
  std::mt19937
      random_number_generator;  // 32-bit unsigned random number generator
};

class RandomInclusiveInt : public RandomGenerator {
 public:
  RandomInclusiveInt() = default;

  RandomInclusiveInt(Seed const &random_seed);

  uint32_t generate(uint32_t min, uint32_t max) override;
};
}  // namespace gram

#endif  // GRAMTOOLS_RANDOM_HPP
