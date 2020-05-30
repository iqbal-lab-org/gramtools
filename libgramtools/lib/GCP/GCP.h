#ifndef GCP_H
#define GCP_H

#include <algorithm>
#include <map>
#include <cmath>
#include <iostream>
#include <vector>
#include <random>

namespace GCP {

using GenotypeConfidence = double;
using GenotypePercentile = double;

/* ________Model________ */

/**
* Class responsible for producing data used by a Genotyper.
* The client is responsible for implementing `produce_data`;
* eg in Pandora, we model k-mer coverage with a negative binomial distribution.
*/
template<typename ModelData>
class Model {
public:
    explicit Model(uint32_t seed = 42) : random_number_generator(seed) {}

    /**
     * To be inherited and implemented by the client.
     * Note that you have access to a random number generator.
     * As the ModelData returned is a template, this class can model anything.
     */
    virtual ModelData produce_data() = 0;

    // destructor
    virtual ~Model() = default;

    // disabling copy and move ctor, and assignment op (always good to be extra safe in C++)
    Model(const Model &other) = delete;

    Model(Model &&other) = delete;

    Model &operator=(const Model &other) = delete;

    Model &operator=(Model &&other) = delete;

protected:
    std::default_random_engine random_number_generator;
};


/* ________Simulator________ */

/**
* Class responsible for simulating genotype confidences given a Model and a Genotyper.
*/
template<typename ModelData, typename Genotyper>
class Simulator {
public:
    // ctor/dtor
    explicit Simulator(Model<ModelData> *model) : model(model) {}

    virtual ~Simulator() = default;

    // main method
    std::vector<GenotypeConfidence> simulate(uint32_t iterations) {
        std::vector<GenotypeConfidence> genotype_confidences;
        genotype_confidences.reserve(iterations);

        for (uint32_t iterations_done = 0; iterations_done < iterations; ++iterations_done) {
            ModelData model_data = model->produce_data();
            Genotyper gtyped(model_data);
            genotype_confidences.push_back(gtyped.get_genotype_confidence());
        }

        std::sort(genotype_confidences.begin(), genotype_confidences.end());
        return genotype_confidences;
    }


    // disabling copy and move ctor, and assignment op (always good to be extra safe in C++)
    Simulator(const Simulator &other) = delete;

    Simulator(Simulator &&other) = delete;

    Simulator &operator=(const Simulator &other) = delete;

    Simulator &operator=(Simulator &&other) = delete;

private:
    Model<ModelData> *model;
};


/* ________Percentiler________ */

class NotEnoughData : public std::runtime_error {
    using std::runtime_error::runtime_error;
};


/**
* Class responsible for assigning confidence percentiles to raw genotype confidences.
* This class is decoupled from the others (to satisfy https://github.com/leoisl/GCP/issues/6)
*/
class Percentiler {
public:
    /**
     * Builds a GenotypeConfidencePercentiler from a vector of simulated genotype confidences
     */
    explicit Percentiler(const std::vector<GenotypeConfidence> &input_entries) {

        bool enough_data = input_entries.size() >= 2;
        if (not enough_data) {
            throw NotEnoughData("Please provide at least two simulated genotype confidences.");
        }

        // Populate map of confidence -> percentile
        auto cur_entry = input_entries.begin();
        while (cur_entry != input_entries.end()){
            // hi is the first entry greater than the cur_entry (or the end iterator)
            auto hi = std::upper_bound(input_entries.begin(), input_entries.end(), *cur_entry);
            auto cur_percentile = iterator_to_percentile(cur_entry, input_entries);
            bool is_single_copy = cur_entry == hi - 1;
            if (is_single_copy) entries[*cur_entry] = cur_percentile;
            else {
                // Case: multiple identical entries, take average
                auto hi_percentile = iterator_to_percentile(hi - 1, input_entries);
                entries[*cur_entry] = cur_percentile + (hi_percentile - cur_percentile) / 2;
            }
            cur_entry = hi;
        }
    }


    /**
     * Get the confidence percentile given a genotype confidence.
     */
    GenotypePercentile get_confidence_percentile(GenotypeConfidence query) {
        auto lower_bound = entries.upper_bound(query);
        if (lower_bound == entries.end()) return 100.0;
        if (lower_bound->first == query) return lower_bound->second;
        if (lower_bound == entries.begin()) return 0.0;

        // Case: need to interpolate between two known confidence/percentile pairs
        auto hi = lower_bound;
        --lower_bound;
        return linear_interpolation(lower_bound->first, hi->first, lower_bound->second,
                hi->second, query);
    }

    // destructor
    virtual ~Percentiler() = default;

    // disabling copy and move ctor, and assignment op (always good to be extra safe in C++)
    Percentiler(const Percentiler &other) = delete;

    Percentiler(Percentiler &&other) = delete;

    Percentiler &operator=(const Percentiler &other) = delete;

    Percentiler &operator=(Percentiler &&other) = delete;

private:
    std::map<GenotypeConfidence, GenotypePercentile> entries;
    using GC_it = std::vector<GenotypeConfidence>::const_iterator;

    GenotypePercentile iterator_to_percentile(GC_it const &it,
            std::vector<GenotypeConfidence> const& input_entries) {
        auto rank = std::distance(input_entries.begin(), it) + 1;
        return 100. * rank / input_entries.size();
    };

    static double linear_interpolation(double x1, double x2, double y1, double y2, double x) {
        auto slope = (y2 - y1) / (x2 - x1);
        return y1 + slope * (x - x1);
    }
};

} // end namespace GCP
#endif // GCP_H
