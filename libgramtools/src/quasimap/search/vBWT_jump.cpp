#include "quasimap/search/vBWT_jump.hpp"

SA_Interval gram::get_allele_marker_sa_interval(const Marker &allele_marker_char,
                                          const PRG_Info &prg_info) {
    const auto alphabet_rank = prg_info.fm_index.char2comp[allele_marker_char];
    const auto start_sa_index = prg_info.fm_index.C[alphabet_rank];

    SA_Index end_sa_index;
    // Case: the current marker is not the last element of the alphabet.
    if (alphabet_rank < prg_info.fm_index.sigma - 1){
        // - 1 because we get positioned at the first position whose suffix starts with the next element in the alphabet
        // and we want to exclude that because we work with inclusive rank queries in backward search.
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
 * Add to the variant path taken, either i)the site and allele marker or ii)the allele marker only
 *
 * Note first that whenever we enter a variant site, we commit a site marker and an (invalid) allele marker to the
 * created search state.
 *
 * So if the variant site path is empty, we have not done this; ie, we started mapping inside the site: case i).
 * If it is not empty, we started mapping outside, so we update the allele id only.
 *
 *
 * @note ABOUT kmer index serialisation
 * Note first that when we enter a var site, we set the `SearchVariantSiteState` enum to 'within_variant'.
 * Now kmer index serialisation is lossy: we lose the `SearchVariantSiteState` if it is set: ie
 * when re-loading the index, it goes to `unknown`.
 * Thus we may have mapped into a site, and reached the kmer end, and serialised that to disk.
 * And thus we cannot use this enum to ascertain if we have previously entered the site.
 */
void update_variant_site_path(SearchState &affected_search_state,
                              const uint64_t allele_id,
                              const Marker site_ID) {
    // Anytime you enter a site, it gets pushed to `traversing_path`
    // If the latter is empty, we have not seen the site entry
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
SearchState exiting_site_search_state(const VariantLocus &locus, const SearchState &current_search_state,
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


SearchStates gram::process_markers_search_states(const SearchStates &old_search_states,
                                                 const PRG_Info &prg_info) {
    SearchStates new_search_states = old_search_states;
    SearchStates all_markers_new_search_states;
    for (const auto &search_state: old_search_states) {
        auto markers_search_states = process_markers_search_state(search_state, prg_info);
        all_markers_new_search_states.splice(all_markers_new_search_states.end(),
                                             markers_search_states);
    }
    new_search_states.splice(new_search_states.end(),
                             all_markers_new_search_states);
    return new_search_states;
}


SearchStates gram::process_markers_search_state(const SearchState &current_search_state,
                                                const PRG_Info &prg_info) {
    // A vector of the `VariantLocus`s that need to be processed
    auto marker_targets = left_markers_search(current_search_state,
                                              prg_info);
    if (marker_targets.empty())
        return SearchStates{};

    target_m const& target_map = prg_info.coverage_graph.target_map;
    SearchStates markers_search_states = {};
    Locus_and_SearchStates extension_targets;
    Locus_and_SearchStates to_process_targets;

    // Add the current search state to the locus for extension of the `SearchState`
    for (auto& marker_target : marker_targets) to_process_targets.push_back({marker_target, current_search_state});

    // In the loop we must respect the following contract:
    // - Each new target has a search state that needs to be committed
    // - A locus is deemed processed, and is thus not processed again, if it is a site exit point
    while (!to_process_targets.empty()) {
        auto const to_process_target = to_process_targets.back();
        to_process_targets.pop_back();
        auto const& target_locus = to_process_target.first;
        auto const& search_state = to_process_target.second;

        // Get the new search state
        if (is_site_marker(target_locus.first)) {
            auto new_target = extend_targets_site_exit(target_locus, search_state, prg_info);
            extension_targets = Locus_and_SearchStates{new_target};
        }

        else {
            extension_targets = extend_targets_site_entry(target_locus, search_state, prg_info);
        }

        // Commit the new target search states
        // and the Loci if they are non site
        for (auto& new_target : extension_targets){
            markers_search_states.push_back(new_target.second);
            auto const& site_ID = new_target.first.first;
            // Add the new target for processing if it is a site entry point (=even marker)
            if (is_allele_marker(site_ID) && site_ID != 0) to_process_targets.push_back(new_target);
        }

    }
    return markers_search_states;
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
        if (is_allele_marker(target_locus.first)){
            if (prg_info.last_allele_positions.at(target_locus.first) !=
                prg_index - 1) target_locus.first--;
        }
        markers_search_results.push_back(target_locus);
    }

    return markers_search_results;
}

Locus_and_SearchState gram::extend_targets_site_exit(VariantLocus const& target_locus, SearchState const& search_state,
                                                     PRG_Info const& prg_info){
    // Make copies that we will modify
    auto new_locus = target_locus;
    auto& site_marker = new_locus.first;
    auto& allele_id = new_locus.second;
    auto new_search_state = search_state;

    auto& target_map = prg_info.coverage_graph.target_map;
    auto prev_id = site_marker;
    while(true){
        // update the SearchState.
        new_search_state = exiting_site_search_state(new_locus, new_search_state, prg_info);

        if (target_map.find(site_marker) == target_map.end()) break; // Nothing to do: no adjacent markers
        // update the VariantLocus.
        auto target_markers = target_map.at(site_marker);
        assert(target_markers.size() == 1); // A site entry point should not point to more than one other marker
        site_marker = target_markers.back().ID;

        if (is_allele_marker(site_marker)) break; // We have an allele marker (site exit); this will get processed separately
        else { // A double exit
            auto parent_site = prg_info.coverage_graph.par_map.at(prev_id);
            // The targeted double exit should be correspondingly well recorded in the parental map
            assert(parent_site.first == site_marker);
            allele_id = parent_site.second;
            prev_id = site_marker;
        }
    }

    return std::make_pair(new_locus, new_search_state);
}

Locus_and_SearchStates gram::extend_targets_site_entry(VariantLocus const& target_locus, SearchState const& search_state,
                                                       PRG_Info const& prg_info){
    // First, simply get a new search state and flag the argument locus as being 'done with'
    auto new_locus = target_locus;
    auto& variant_marker = new_locus.first;
    auto& allele_id = new_locus.second;

    auto original_allele_marker = variant_marker;
    variant_marker = 1; // A dummy value guaranteeing the target will not be processed again in the calling loop.
    Locus_and_SearchStates extensions;
    auto new_search_state = entering_site_search_state(target_locus.first, search_state, prg_info);
    extensions.push_back({new_locus, new_search_state});

    // Now look for extensions
    auto& target_map = prg_info.coverage_graph.target_map;
    bool has_targets = target_map.find(original_allele_marker) != target_map.end();

    if (!has_targets) return extensions;

    // Traverse each target in the map and add it as an extension
    for(auto& mapped_target : target_map.at(original_allele_marker)){
        if (is_site_marker(mapped_target.ID)){ // Case: direct deletion
            assert(mapped_target.direct_deletion_allele != 0);
            VariantLocus site_exit_locus{mapped_target.ID, mapped_target.direct_deletion_allele};
            // Add a fake traversing element that gets specified in the call to `extend_targets_site_exit`
            auto direct_deletion_search_state = new_search_state;
            direct_deletion_search_state.traversing_path.push_back(VariantLocus{mapped_target.ID, ALLELE_UNKNOWN});
            auto extended_target = extend_targets_site_exit(site_exit_locus, new_search_state, prg_info);
            extensions.push_back(extended_target);
        }

        else { // Case : double entry
            VariantLocus site_entry_locus{mapped_target.ID, ALLELE_UNKNOWN};
            extensions.push_back({site_entry_locus, new_search_state});
        }
    }
    return extensions;
}

