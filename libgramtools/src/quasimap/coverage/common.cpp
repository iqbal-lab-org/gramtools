#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/nondet_random.hpp>

#include "quasimap/coverage/allele_sum.hpp"
#include "quasimap/coverage/allele_base.hpp"
#include "quasimap/coverage/grouped_allele_counts.hpp"

#include "quasimap/coverage/common.hpp"


SA_Interval random_select_sa_interval(const SearchStates &search_states,
                                      const uint32_t &random_seed) {
    uint32_t actual_seed = 0;
    if (random_seed == 0) {
        boost::random_device seed_generator;
        actual_seed = seed_generator();
    } else {
        actual_seed = random_seed;
    }
    boost::mt19937 random_number_generator(actual_seed);
    boost::uniform_int<> range(0, (int) search_states.size() - 1);

    using Generator = boost::variate_generator<boost::mt19937, boost::uniform_int<>>;
    Generator generate_random_number(random_number_generator, range);

    auto it_advance_steps = generate_random_number();
    auto it = search_states.begin();
    std::advance(it, it_advance_steps);
    const auto &search_state = *it;
    return search_state.sa_interval;
}


SearchStates filter_for_sa_interval(const SA_Interval &target_sa_interval,
                                    const SearchStates &search_states) {
    SearchStates new_search_states = {};
    for (const auto &search_state: search_states) {
        if (search_state.sa_interval != target_sa_interval)
            continue;
        new_search_states.push_back(search_state);
    }
    return new_search_states;
}


void coverage::record::search_states(Coverage &coverage,
                                     const SearchStates &search_states,
                                     const uint64_t &read_length,
                                     const PRG_Info &prg_info,
                                     const uint32_t &random_seed) {
    auto sa_interval = random_select_sa_interval(search_states, random_seed);
    auto filtered_search_states = filter_for_sa_interval(sa_interval, search_states);

    coverage::record::allele_sum(coverage, filtered_search_states);
    coverage::record::grouped_allele_counts(coverage, filtered_search_states);
    coverage::record::allele_base(coverage, filtered_search_states, read_length, prg_info);
}


void coverage::dump::all(const Coverage &coverage,
                         const Parameters &parameters) {
    coverage::dump::allele_sum(coverage, parameters);
    coverage::dump::allele_base(coverage, parameters);
    coverage::dump::grouped_allele_counts(coverage, parameters);
}


Coverage coverage::generate::empty_structure(const PRG_Info &prg_info) {
    Coverage coverage = {};
    coverage.allele_sum_coverage = coverage::generate::allele_sum_structure(prg_info);
    coverage.allele_base_coverage = coverage::generate::allele_base_structure(prg_info);
    coverage.grouped_allele_counts = coverage::generate::grouped_allele_counts(prg_info);
    return coverage;
}