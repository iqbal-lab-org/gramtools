#include "common/random.hpp"

namespace gram {
RandomInclusiveInt::RandomInclusiveInt(Seed const &random_seed) {
  SeedSize master_seed;
  if (random_seed.has_value())
    master_seed = random_seed.value();
  else {
    std::random_device seed_generator;
    master_seed = seed_generator();
  }
  this->random_number_generator.seed(master_seed);
}

uint32_t RandomInclusiveInt::generate(uint32_t min, uint32_t max) {
  std::uniform_int_distribution<uint32_t> range(min, max);
  return range(random_number_generator);
}
}  // namespace gram
