#ifndef GRAMTOOLS_RANDOM_HPP
#define GRAMTOOLS_RANDOM_HPP

#include <cstdint>

namespace gram {

    /**
     * Abstract base class used for mocking in unit tests
     */
    class RandomGenerator {
    public:
        virtual ~RandomGenerator() {};

        virtual uint32_t generate(uint32_t min, uint32_t max) const = 0;
    };

    class RandomInclusiveInt : public RandomGenerator {
    public:
        RandomInclusiveInt() = default;

        RandomInclusiveInt(uint32_t const &random_seed);

        uint32_t generate(uint32_t min, uint32_t max) const override;

    private:
        uint32_t random_seed;
    };
}

#endif //GRAMTOOLS_RANDOM_HPP
