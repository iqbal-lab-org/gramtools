#include "quasimap/search/vBWT_jump.hpp"

SA_Interval gram::get_allele_marker_sa_interval(const Marker &allele_marker_char,
                                                const PRG_Info &prg_info) {
    const auto alphabet_rank = prg_info.fm_index.char2comp[allele_marker_char];
    const auto start_sa_index = prg_info.fm_index.C[alphabet_rank];

    SA_Index end_sa_index;
    // Case: the current marker is not the last element of the alphabet.
    if (alphabet_rank < prg_info.fm_index.sigma - 1) {
        // Below: use - 1 because we get positioned at the first position whose suffix starts with the next element in the alphabet
        // note: the rank query itself is exclusive, so at backward search time this will get +1 again.
        end_sa_index = prg_info.fm_index.C[alphabet_rank + 1] - 1;
    }
        // Case: it is the last element of the alphabet
    else end_sa_index = prg_info.fm_index.size() - 1;

    return SA_Interval{start_sa_index, end_sa_index};
}


/**
 * Given an allele (=even) marker SA interval, make one search state containing them all.
 * We will record the allele id later, if and when each allele gets out.
 * @param current_search_state: gets passed so as to chain onto list of previously traversed VariantLoci
 */
SearchState get_allele_search_state(const Marker &allele_marker,
                                    const SA_Interval &allele_marker_sa_interval,
                                    const SearchState &current_search_state,
                                    const PRG_Info &prg_info) {


    auto site_boundary_marker = allele_marker - 1;
    SearchState search_state = current_search_state;
    search_state.sa_interval.first = allele_marker_sa_interval.first;
    search_state.sa_interval.second = allele_marker_sa_interval.second;

    search_state.variant_site_state
            = SearchVariantSiteState::within_variant_site;

    search_state.traversing_path.push_back(VariantLocus{site_boundary_marker, ALLELE_UNKNOWN});

    return search_state;
}

/**
 * Deals with a read mapping into a variant site's end point.
 * The SA index of each allele's end gets added as a new SearchState.
 * Because a variant site end is found, the read needs to be able to map through all alleles of this site.
 */
SearchState entering_site_search_state(const Marker &allele_marker,
                                       const SearchState &current_search_state,
                                       const PRG_Info &prg_info) {

    // Get full SA interval of the corresponding allele marker.
    auto allele_marker_sa_interval =
            get_allele_marker_sa_interval(allele_marker,
                                          prg_info);
    // Get one `SearchState` for the whole site, with populated cache.
    auto allele_search_state = get_allele_search_state(allele_marker,
                                                       allele_marker_sa_interval,
                                                       current_search_state,
                                                       prg_info);

    return allele_search_state;
}

/**
 * Add to the variant path taken the site and allele markers
 * They can be present at the back of the `traversing_path` in which case we appropriately set the allele_id
 * and move the marker to the `traversed_path`
 *
 * @note ABOUT kmer index serialisation
 * When we enter a var site, we set the `SearchVariantSiteState` enum to 'within_variant'.
 * Kmer index serialisation is lossy: we lose the `SearchVariantSiteState` if it is set: ie
 * when re-loading the index, it goes to `unknown`.
 * Thus we should not rely on this enum to * ascertain if we have previously entered the site.
 */
void update_variant_site_path(SearchState &affected_search_state,
                              const uint64_t allele_id,
                              const Marker site_ID) {
    // Anytime you enter a site, it gets pushed to `traversing_path`
    // If the latter is empty, we have not seen the site entry (ie, we started mapping inside site)
    bool started_in_site = affected_search_state.traversing_path.empty();
    if (started_in_site) { // Case: make new site/allele pair
        affected_search_state.traversed_path.push_back(VariantLocus{site_ID, allele_id});
    } else { // Case: add the allele id to existing site/allele pair
        auto existing_locus = affected_search_state.traversing_path.back();
        // Make sure we're recording leaving the right site
        assert(existing_locus.first == site_ID);
        assert(existing_locus.second == ALLELE_UNKNOWN);
        existing_locus.second = allele_id;
        affected_search_state.traversed_path.push_back(existing_locus);
        affected_search_state.traversing_path.pop_back();
    }
}

/**
 * Deals with a read mapping leaving a variant site.
 * Create a new `SearchState` with SA interval the index of the site variant's entry point.
 */
SearchState exiting_site_search_state(const VariantLocus &locus,
                                      const SearchState &current_search_state,
                                      const PRG_Info &prg_info) {

    SearchState new_search_state = current_search_state;

    Marker site_marker{locus.first};
    AlleleId allele_id{locus.second};

    update_variant_site_path(new_search_state, allele_id, site_marker);

    auto alphabet_rank = prg_info.fm_index.char2comp[site_marker];
    SA_Index site_index = prg_info.fm_index.C[alphabet_rank];

    new_search_state.sa_interval = SA_Interval{site_index, site_index};
    new_search_state.variant_site_state
            = SearchVariantSiteState::outside_variant_site;

    return new_search_state;
}

MarkersSearchResults gram::left_markers_search(const SearchState &search_state,
                                               const PRG_Info &prg_info) {
    MarkersSearchResults markers_search_results;

    const auto &sa_interval = search_state.sa_interval;

    for (int index = sa_interval.first; index <= sa_interval.second; index++) {
        if (prg_info.bwt_markers_mask[index] == 0) continue;

        auto prg_index = prg_info.fm_index[index];
        VariantLocus target_locus = prg_info.coverage_graph.random_access[prg_index].target;
        // Convert the target to a site ID if it is an allele ID that points to the beginning of the site
        // (ie, it is not the last allele)
        if (is_allele_marker(target_locus.first)) {
            if (prg_info.last_allele_positions.at(target_locus.first) !=
                prg_index - 1)
                target_locus.first--;
        }
        markers_search_results.push_back(target_locus);
    }

    return markers_search_results;
}

SearchStates gram::process_markers_search_states(const SearchStates &old_search_states,
                                                 const PRG_Info &prg_info) {
    SearchStates new_search_states = old_search_states;
    SearchStates all_markers_new_search_states;
    for (const auto &search_state: old_search_states) {
        auto markers_search_states = search_state_vBWT_jumps(search_state, prg_info);
        all_markers_new_search_states.splice(all_markers_new_search_states.end(),
                                             markers_search_states);
    }
    new_search_states.splice(new_search_states.end(),
                             all_markers_new_search_states);
    return new_search_states;
}


SearchStates gram::search_state_vBWT_jumps(const SearchState &current_search_state,
                                           const PRG_Info &prg_info) {
    // A vector of the `VariantLocus`s that need to be processed
    auto marker_targets = left_markers_search(current_search_state,
                                              prg_info);
    if (marker_targets.empty())
        return SearchStates{};

    target_m const &target_map = prg_info.coverage_graph.target_map;
    SearchStates markers_search_states = {};
    Locus_and_SearchStates extension_targets;
    Locus_and_SearchStates to_process_targets;

    // Add the current search state to each locus; each will be extended independently
    for (auto &marker_target : marker_targets) to_process_targets.push_back({marker_target, current_search_state});

    // In the loop we must respect the following contract:
    // - Each new target has a search state that says if it needs to be committed
    // - A locus is deemed processed, and is thus not processed again, if it is a site exit point
    while (!to_process_targets.empty()) {
        auto const to_process_target = to_process_targets.back();
        to_process_targets.pop_back();
        auto const &target_locus = to_process_target.locus;
        auto const &search_state = to_process_target.search_state;

        // Get the new targets
        if (is_site_marker(target_locus.first)) {
            auto new_target = extend_targets_site_exit(target_locus, search_state, prg_info);
            extension_targets = Locus_and_SearchStates{new_target};
        } else {
            extension_targets = extend_targets_site_entry(target_locus, search_state, prg_info);
        }

        // Commit the new target search states and loci
        for (auto &new_target : extension_targets) {
            // Does the target need to be backward searched?
            if (new_target.commit_me) markers_search_states.push_back(new_target.search_state);

            // Does the target need to be processed further due to adjacent variant markers?
            auto const &site_ID = new_target.locus.first;
            if (site_ID != 0) to_process_targets.push_back(new_target);
        }
    }
    return markers_search_states;
}



Locus_and_SearchState gram::extend_targets_site_exit(VariantLocus const &target_locus, SearchState const &search_state,
                                                     PRG_Info const &prg_info) {

    VariantLocus next_target = target_locus;
    auto site_marker = next_target.first;
    bool commit_me{true};
    auto &target_map = prg_info.coverage_graph.target_map;

    // update the SearchState.
    auto new_search_state = exiting_site_search_state(target_locus, search_state, prg_info);
    // Signal we do not want to process the locus further, by default.
    next_target = VariantLocus{0,0};

    while (target_map.find(site_marker) != target_map.end()) {

        auto target_markers = target_map.at(site_marker);
        assert(target_markers.size() == 1); // A site entry point should not point to more than one other marker

        // Make new_locus point to the next site marker that needs a vBWT jump
        auto next_site_marker = target_markers.back().ID;

        if (is_allele_marker(next_site_marker)){ // An exit followed by an entry
            next_target = VariantLocus{next_site_marker, 0};
            commit_me = false; // There will be no sequence to extend into, so we will not extend this SearchState
            break;  // The site entry will get processed separately by the calling function
        }
        else { // A double exit
            // Sanity check: the targeted double exit should be correspondingly well recorded in the parental map
            auto parent_site = prg_info.coverage_graph.par_map.at(site_marker);
            assert(parent_site.first == next_site_marker);

            // update the SearchState.
            auto allele_id = parent_site.second;
            new_search_state = exiting_site_search_state(VariantLocus{next_site_marker, allele_id},
                                                            new_search_state, prg_info);
            site_marker = next_site_marker;
        }
    }
    return Locus_and_SearchState{next_target, new_search_state, commit_me};
}

Locus_and_SearchStates
gram::extend_targets_site_entry(VariantLocus const &target_locus, SearchState const &search_state,
                                PRG_Info const &prg_info) {
    Locus_and_SearchStates extensions;
    VariantLocus next_target;

    auto variant_marker = target_locus.first;

    // First, simply ready the search state for mapping into the site, and flag the locus as being 'done with'
    auto new_search_state = entering_site_search_state(target_locus.first, search_state, prg_info);
    next_target = VariantLocus{0,0};
    extensions.push_back({next_target, new_search_state, true});

    // Now look for extensions
    auto &target_map = prg_info.coverage_graph.target_map;
    bool has_targets = target_map.find(variant_marker) != target_map.end();

    if (!has_targets) return extensions;

    // Traverse each target in the map and add it as an extension
    for (auto &mapped_target : target_map.at(variant_marker)) {
        if (is_site_marker(mapped_target.ID)) { // Case: direct deletion
            assert(mapped_target.direct_deletion_allele != 0);
            VariantLocus site_exit_locus{mapped_target.ID, mapped_target.direct_deletion_allele};

            // Add a traversing element for the direct deletion whose allele ID will later get specified,
            // Even though we know the allele ID right here.
            // This maintains compatibility with mapping occurring without adjacent markers.
            auto direct_deletion_search_state = new_search_state;
            direct_deletion_search_state.traversing_path.push_back(VariantLocus{mapped_target.ID, ALLELE_UNKNOWN});

            extensions.push_back({site_exit_locus, direct_deletion_search_state, false});

        } else { // Case : double entry
            VariantLocus site_entry_locus{mapped_target.ID, ALLELE_UNKNOWN};
            extensions.push_back({site_entry_locus, new_search_state, false});
        }
    }
    return extensions;
}

