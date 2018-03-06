#include <unordered_set>

#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/nondet_random.hpp>

#include "quasimap/coverage/allele_sum.hpp"
#include "quasimap/coverage/allele_base.hpp"
#include "quasimap/coverage/grouped_allele_counts.hpp"

#include "quasimap/coverage/common.hpp"


SearchState random_select_search_state(const SearchStates &search_states,
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
    const auto &choice = *it;
    return choice;
}


bool check_allele_encapsulated(const SearchState &search_state,
                               const uint64_t &read_length,
                               const PRG_Info &prg_info) {
    bool single_allele_path = search_state.variant_site_path.size() == 1;
    bool start_within_allele = search_state.variant_site_state == SearchVariantSiteState::within_variant_site;

    bool end_within_allele = true;
    for (SA_Index sa_index = search_state.sa_interval.first;
         sa_index <= search_state.sa_interval.second;
         ++sa_index) {
        auto start_index = prg_info.fm_index[sa_index];
        auto start_site_marker = prg_info.sites_mask[start_index];
        auto start_allele_id = prg_info.allele_mask[start_index];

        auto end_index = start_index + read_length - 1;
        assert(end_index < prg_info.encoded_prg.size());

        auto end_site_marker = prg_info.sites_mask[end_index];
        auto end_allele_id = prg_info.allele_mask[end_index];
        end_within_allele = (start_site_marker == end_site_marker)
                            and (start_allele_id == end_allele_id);
        if (not end_within_allele)
            break;
    }

    return single_allele_path and start_within_allele and end_within_allele;
}


bool multiple_allele_encapsulated(const SearchState &search_state,
                                  const uint64_t &read_length,
                                  const PRG_Info &prg_info) {
    bool allele_encapsulated = check_allele_encapsulated(search_state,
                                                         read_length,
                                                         prg_info);
    bool multiple_mappings = search_state.sa_interval.second - search_state.sa_interval.first > 0;
    return allele_encapsulated and multiple_mappings;
}


SearchState random_select_single_mapping(const SearchState &search_state,
                                         const uint32_t &random_seed) {
    uint32_t actual_seed = 0;
    if (random_seed == 0) {
        boost::random_device seed_generator;
        actual_seed = seed_generator();
    } else {
        actual_seed = random_seed;
    }
    boost::mt19937 random_number_generator(actual_seed);

    boost::uniform_int<> range((int) search_state.sa_interval.first, (int) search_state.sa_interval.second);
    using Generator = boost::variate_generator<boost::mt19937, boost::uniform_int<>>;
    Generator generate_random_number(random_number_generator, range);
    auto selected_sa_index = (SA_Index) generate_random_number();

    SearchState single_mapping_search_state = search_state;
    single_mapping_search_state.sa_interval.first = selected_sa_index;
    single_mapping_search_state.sa_interval.second = selected_sa_index;
    return single_mapping_search_state;
}


void coverage::record::search_states(Coverage &coverage,
                                     const SearchStates &search_states,
                                     const uint64_t &read_length,
                                     const PRG_Info &prg_info,
                                     const uint32_t &random_seed) {
    auto search_state = random_select_search_state(search_states, random_seed);
    bool has_path = not search_state.variant_site_path.empty();
    if (not has_path)
        return;

    if (multiple_allele_encapsulated(search_state, read_length, prg_info)) {
        search_state = random_select_single_mapping(search_state,
                                                    random_seed);
    }
    SearchStates chosen_search_states = {search_state};

    coverage::record::allele_sum(coverage, chosen_search_states);
    coverage::record::grouped_allele_counts(coverage, chosen_search_states);
    coverage::record::allele_base(coverage, chosen_search_states, read_length, prg_info);
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