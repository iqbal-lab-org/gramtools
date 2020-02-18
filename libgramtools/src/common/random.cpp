#include "common/random.hpp"
#include <boost/random.hpp>
#include <boost/nondet_random.hpp>

namespace gram {
    RandomInclusiveInt::RandomInclusiveInt(uint32_t const &random_seed) {
        if (random_seed == 0) {
            boost::random_device seed_generator;
            this->random_seed = seed_generator();
        } else this->random_seed = random_seed;
    }


    uint32_t RandomInclusiveInt::generate(uint32_t min, uint32_t max) const {
        boost::mt19937 random_number_generator(random_seed);

        boost::uniform_int<> range(min, max);
        using Generator = boost::variate_generator<boost::mt19937, boost::uniform_int<>>;
        Generator generate_random_number(random_number_generator, range);
        return generate_random_number();
    };
}
