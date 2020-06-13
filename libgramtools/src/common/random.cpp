#include "common/random.hpp"
#include <random>

namespace gram {
RandomInclusiveInt::RandomInclusiveInt(uint32_t const &random_seed) {
  if (random_seed == 0) {
    std::random_device seed_generator;
    this->random_seed = seed_generator();
  } else
    this->random_seed = random_seed;
}

uint32_t RandomInclusiveInt::generate(uint32_t min, uint32_t max) const {
  std::mt19937_64 random_number_generator(random_seed);

  std::uniform_int_distribution<std::mt19937_64::result_type> range(min, max);
  return range(random_number_generator);
}
}  // namespace gram
