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

LocusFinder::LocusFinder(SearchState const search_state, info_ptr prg_info)
: search_state(search_state), prg_info(prg_info){
    assign_traversing_loci();
    assign_traversed_loci();
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
      auto allele_id = node_access.node->get_allele();

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

void MappingInstanceSelector::add_searchstate(SearchState const& ss){
    LocusFinder l{ss, prg_info};
    coverage_struct c{ss, l.unique_loci};
    usps[l.base_sites].push_back(c);
}


bool gram::check_allele_encapsulated(const SearchState &search_state,
                                     const uint64_t &read_length,
                                     const PRG_Info &prg_info) {
    bool single_allele_path = search_state.traversed_path.size() == 1;
    bool start_within_allele = search_state.variant_site_state
                               == SearchVariantSiteState::within_variant_site;

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


bool gram::multiple_allele_encapsulated(const SearchState &search_state,
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


uint64_t gram::random_int_inclusive(const uint64_t &min,
                                    const uint64_t &max,
                                    const uint64_t &random_seed) {
    uint64_t actual_seed = 0;
    if (random_seed == 0) {
        boost::random_device seed_generator;
        actual_seed = seed_generator();
    } else {
        actual_seed = random_seed;
    }
    boost::mt19937 random_number_generator((uint32_t) actual_seed);

    boost::uniform_int<> range((uint32_t) min, (uint32_t) max);
    using Generator = boost::variate_generator<boost::mt19937, boost::uniform_int<>>;
    Generator generate_random_number(random_number_generator, range);
    return (uint64_t) generate_random_number();
}


uint32_t gram::count_nonvariant_search_states(const SearchStates &search_states) {
    uint32_t count = 0;
    for (const auto &search_state: search_states) {
        bool has_path =  search_state_has_path(search_state);
        if (not has_path)
            ++count;
    }
    return count;
}

void add_path_sites(SitePath& site_path, VariantSitePath const& copy_from){
    for (auto const& entry: copy_from){
        Marker site = entry.first;
        if (site_path.find(site) != site_path.end()){
            throw std::logic_error(
                     "Trying to add site marker to site path but it was there already.\n "
                     "Site is marked as being traversed twice by the searchstate."
            );
        }
        site_path.insert(site);
    }
}

SitePath gram::get_path_sites(const SearchState &search_state) {
    SitePath site_path = {};
    add_path_sites(site_path, search_state.traversing_path);
    add_path_sites(site_path, search_state.traversed_path);
    return site_path;
}


uniqueSitePaths gram::get_unique_site_paths(const SearchStates &search_states) {
    uniqueSitePaths all_path_sites;
    for (const auto &search_state: search_states) {
        bool has_path =  search_state_has_path(search_state);
        if (not has_path)
            continue;

        auto site_path = get_path_sites(search_state);
        all_path_sites[site_path].push_back(search_state);
    }
    return all_path_sites;
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
 *  * The read maps through a single site/allele combination, multiple times. This means the read is **encapsulated** inside
 *  two different alleles of the same site. Randomly select and return only one of those.
 */
SearchStates selection(const SearchStates &search_states,
                       const uint64_t &read_length,
                       const PRG_Info &prg_info,
                       const uint32_t &random_seed) {
    uint64_t nonvariant_count = count_nonvariant_search_states(search_states);
    // Extract the unique path sites: vectors of site IDs traversed.
    auto path_sites = get_unique_site_paths(search_states);

    uint64_t count_total_options = nonvariant_count + path_sites.size();
    if (count_total_options == 0)
        return SearchStates{};
    uint64_t selected_option = random_int_inclusive(1, count_total_options, random_seed);

    // If we select a non-variant path, return empty `SearchStates`, leading to no coverage information.
    bool selected_no_path = selected_option <= nonvariant_count;
    if (selected_no_path)
        return SearchStates{};

    // We have selected a variant path.
    uint64_t paths_sites_offset = selected_option - nonvariant_count - 1;
    auto it = path_sites.begin();
    std::advance(it, paths_sites_offset);
    auto selected_path_sites = it->first; // a single vector of site IDs (sites path)

    // Filter all given search states for the single randomly selected vector of site_id (site path).
    // Allele IDs are not considered at this point.
    // Mapping instance interactions within sites are considered later (in grouped allele coverage).
    auto selected_search_states = it->second;
    if (selected_search_states.size() > 1)
        return selected_search_states;

    auto search_state = selected_search_states.front();
    if (multiple_allele_encapsulated(search_state, read_length, prg_info)) {
        search_state = random_select_single_mapping(search_state,
                                                    random_seed);
    }
    return SearchStates{search_state};
}


void coverage::record::search_states(Coverage &coverage,
                                     const SearchStates &search_states,
                                     const uint64_t &read_length,
                                     const PRG_Info &prg_info,
                                     const uint32_t &random_seed) {
    SearchStates selected_search_states = selection(search_states,
                                                    read_length,
                                                    prg_info,
                                                    random_seed);

    coverage::record::allele_sum(coverage, selected_search_states);
    coverage::record::grouped_allele_counts(coverage, selected_search_states);
    coverage::record::allele_base(coverage, selected_search_states, read_length, prg_info);
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