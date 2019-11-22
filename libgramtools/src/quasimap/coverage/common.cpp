#include <boost/random.hpp>
#include <boost/nondet_random.hpp>

#include "quasimap/coverage/allele_sum.hpp"
#include "quasimap/coverage/allele_base.hpp"
#include "quasimap/coverage/grouped_allele_counts.hpp"

#include "quasimap/coverage/common.hpp"


using namespace gram;

/**
 * A `SearchState` has a path if either its traversed_path or its traversing_path is non empty
 */
bool search_state_has_path(SearchState const& search_state){
    bool has_path = not search_state.traversed_path.empty();
    has_path = has_path || not search_state.traversing_path.empty();
    return has_path;
}

RandomInclusiveInt::RandomInclusiveInt(uint32_t const& random_seed){
    if (random_seed == 0){
        boost::random_device seed_generator;
        this->random_seed = seed_generator();
    }
    else this->random_seed = random_seed;
}


uint32_t RandomInclusiveInt::generate(uint32_t min, uint32_t max){
    boost::mt19937 random_number_generator(random_seed);

    boost::uniform_int<> range(min, max);
    using Generator = boost::variate_generator<boost::mt19937, boost::uniform_int<>>;
    Generator generate_random_number(random_number_generator, range);
    return generate_random_number();
};

LocusFinder::LocusFinder(SearchState const search_state, info_ptr prg_info)
: search_state(search_state), prg_info(prg_info){
    check_site_uniqueness();
    assign_traversing_loci();
    assign_traversed_loci();
}

void LocusFinder::check_site_uniqueness(SearchState const& search_state) {
    auto all_loci = search_state.traversed_path;
    all_loci.insert(all_loci.end(), search_state.traversing_path.begin(), search_state.traversing_path.end());
    SitePath unique_sites;
    for (auto const &entry: all_loci) {
        Marker site = entry.first;
        if (unique_sites.find(site) != unique_sites.end()) {
            throw std::logic_error(
                    "ERROR: A site cannot have been traversed more than once by a read, but this one is marked as such.\n"
            );
        }
        unique_sites.insert(site);
    }
    return;
}

void LocusFinder::assign_nested_locus(VariantLocus const& var_loc, info_ptr info_ptr){
    auto& par_map = info_ptr->coverage_graph.par_map;
    VariantLocus cur_locus = var_loc;
    Marker& cur_marker = cur_locus.first;
   while(true){
      if (used_sites.find(cur_marker) != used_sites.end()) break;
      used_sites.insert(cur_marker);
      unique_loci.insert(cur_locus);

      if (par_map.find(cur_marker) == par_map.end()){
          base_sites.insert(cur_marker); // Add non-nested site marker
          break;
      }
      cur_locus = par_map.at(cur_marker);
   }
   return;
}

void LocusFinder::assign_traversing_loci(SearchState const& search_state, info_ptr prg_info){
    if (search_state.traversing_path.size() == 0) return;
    auto r = search_state.traversing_path.rbegin();
    Marker parent_seed = r->first;
    assert(r->second == ALLELE_UNKNOWN);

    VariantLocus new_locus;
    // Assign the currently traversed alleles
    for (int i = search_state.sa_interval.first; i <= search_state.sa_interval.second; ++i){
      auto prg_pos =  prg_info->fm_index[i];
      auto& node_access = prg_info->coverage_graph.random_access[prg_pos];
      auto allele_id = node_access.node->get_allele_ID();

      new_locus = VariantLocus{parent_seed, allele_id};
      unique_loci.insert(new_locus);
    }

    assign_nested_locus(new_locus, prg_info);

    // TODO: add a check that entries in parent map correspond to entrie in traversin_locus vector
    // auto r = search_state.traversing_path.rbegin();
}

void LocusFinder::assign_traversed_loci(SearchState const& search_state, info_ptr prg_info){
    for (auto const & var_locus : search_state.traversed_path){
        assign_nested_locus(var_locus, prg_info);
    }
}

MappingInstanceSelector::MappingInstanceSelector(SearchStates const search_states, info_ptr prg_info, rand_ptr rand_generator)
    : input_search_states(search_states), usps(), prg_info(prg_info), rand_generator(rand_generator)
{
    add_searchstates();
    int32_t selected_index = random_select_entry();
    if (selected_index >= 0) apply_selection(selected_index);
};

int32_t MappingInstanceSelector::random_select_entry() {
    if (usps.size() == 0) return -1;
    uint32_t nonvariant_count = count_nonvar_search_states();
    uint32_t count_total_options = nonvariant_count + usps.size();

    uint32_t selected_option = rand_generator->generate(1, count_total_options);
    // If we select a non-variant path, no coverage information will get recorded.
    bool no_variants = selected_option <= nonvariant_count;
    if (no_variants) return -1;

    return selected_option - nonvariant_count - 1;  //return 0-based index in map
}

void MappingInstanceSelector::apply_selection(int32_t selected_index){
    auto it = usps.begin();
    std::advance(it, selected_index);
    auto chosen_traversal = it->second;
    selected = SelectedMapping{chosen_traversal.first, chosen_traversal.second};
}

void MappingInstanceSelector::add_searchstate(SearchState const& ss){
    LocusFinder l{ss, prg_info};
    // Create or retrieve the coverage information
    auto& cov_info = usps[l.base_sites];
    // Merge each locus to the existing coverage information
    for (auto& locus : l.unique_loci) cov_info.second.insert(locus);
    cov_info.first.push_back(ss);
}

void MappingInstanceSelector::add_searchstates(SearchStates const& all_ss){
    for (auto const& ss : all_ss){
        if (search_state_has_path(ss)) add_searchstate(ss);
    }
}

uint32_t MappingInstanceSelector::count_nonvar_search_states(SearchStates const& search_states) {
    uint32_t count = 0;
    for (const auto &search_state: search_states) {
        bool has_path =  search_state_has_path(search_state);
        if (not has_path)
            // Add all distinct mappings
            count += (search_state.sa_interval.second - search_state.sa_interval.first + 1);
    }
    return count;
}

/**
 * Random uniform selection of one mapped path in the prg.
 * What gets selected is either:
 * * A single non-variant path through the prg, chosen from all distinct mappings of the read in the prg going through non-variant sites.
 * * A single unique variant path through the prg. Unique is defined as a set of site IDs, irrespective of allele IDs inside sites.
 *
 *  In the second case, there are several sub-cases when **including allele ID traversal**:
 *  * A single site/allele combination path has been selected. Return this.
 *  * The read maps several times through the same set of sites, going through different alleles. Return all of them. Their coverages will all get recorded.
 *  Specifically this will also be used for reporting allele group counts.
 *  * [TODO] The read has horizontal uncertainty: for eg maps inside one site/allele combination twice.
 *  Probably need to randomly select only one mapping instance from those subcases.
 */
SelectedMapping selection(const SearchStates &search_states,
                       const uint64_t &read_length,
                       const PRG_Info &prg_info,
                       const uint32_t &random_seed) {
    RandomInclusiveInt rng{random_seed};
    MappingInstanceSelector m{search_states, &prg_info, &rng};

    // This contains empty containers if we selected a mapping instance in an invariant part of the PRG
    SelectedMapping selected = m.get_selection();

    return selected;
}

void gram::set_allele_ids(SearchStates &search_states,
                          const PRG_Info &prg_info) {

    for (auto &search_state: search_states) {
        // If the variant site path is empty, we cannot have an unknown allele id.
        if (!search_state.traversing_path.empty()) {

            auto last = search_state.traversing_path.back();
            search_state.traversing_path.pop_back();
            assert(search_state.traversing_path.size() == 0);

            for (int SA_pos = search_state.sa_interval.first; SA_pos <= search_state.sa_interval.second; SA_pos++) {
                auto new_search_state = search_state;
                // Find out allele id
                auto text_pos = prg_info.fm_index[SA_pos];
                auto allele_id = prg_info.allele_mask[text_pos];
                last.second = allele_id;

                new_search_state.traversed_path.emplace_back(last);
                new_search_state.sa_interval = SA_Interval{SA_pos, SA_pos};

                if (SA_pos != search_state.sa_interval.second){ // Case: add to the set of search states
                    search_states.emplace_back(new_search_state);
                }
                else search_state = new_search_state; // Case: end of the iteration; modify the search state in place.
            }
        }
    }
}


void coverage::record::search_states(Coverage &coverage,
                                     const SearchStates &search_states,
                                     const uint64_t &read_length,
                                     const PRG_Info &prg_info,
                                     const uint32_t &random_seed) {
    SelectedMapping selected_search_states = selection(search_states,
                                                    read_length,
                                                    prg_info,
                                                    random_seed);

    // If we selected a mapping instance that does not overlap any variant site, there is no coverage to record.
    if (selected_search_states.navigational_search_states.size() == 0) return;

    set_allele_ids(selected_search_states.navigational_search_states, prg_info);

    // [TODO :remove] trick to only run these on non-nested PRG strings.
    if (prg_info.coverage_graph.par_map.size() == 0){
        coverage::record::allele_base(coverage, selected_search_states.navigational_search_states, read_length, prg_info);
    }
    coverage::record::allele_sum(coverage, selected_search_states.equivalence_class_loci);
    coverage::record::grouped_allele_counts(coverage, selected_search_states.equivalence_class_loci);
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