#include <sdsl/suffix_arrays.hpp>
#include "search/search.hpp"


using namespace gram;

/**
 * A caching object used to temporarily store a single search state
 * @see handle_allele_encapsulated_state()
 */
class SearchStateCache {
public:
    SearchState search_state = {};
    bool empty = true;

    void set(const SearchState &search_state) {
        this->search_state = search_state;
        this->empty = false;
    }

    void flush(SearchStates &search_states) {
        if (this->empty)
            return;
        search_states.emplace_back(this->search_state);
        this->empty = true;
    }

    void update_sa_interval_max(const SA_Index &new_sa_interval_max) {
        assert(not this->empty);
        assert(this->search_state.sa_interval.second + 1 == new_sa_interval_max);
        this->search_state.sa_interval.second = new_sa_interval_max;
    }
};


SearchStates gram::handle_allele_encapsulated_state(const SearchState &search_state,
                                                    const PRG_Info &prg_info) {
    bool has_path = not search_state.traversed_path.empty();
    assert(not has_path);

    SearchStates new_search_states = {};
    SearchStateCache cache;

    for (uint64_t sa_index = search_state.sa_interval.first;
         sa_index <= search_state.sa_interval.second;
         ++sa_index) {

        auto prg_index = prg_info.fm_index[sa_index];
        auto site_marker = prg_info.sites_mask[prg_index];
        auto allele_id = prg_info.allele_mask[prg_index];

        bool within_site = site_marker != 0;
        if (not within_site) {
            cache.flush(new_search_states);
            cache.set(SearchState{
                    SA_Interval{sa_index, sa_index},
                    VariantSitePath{},
                    VariantSitePath{},
                    SearchVariantSiteState::outside_variant_site
            });
            cache.flush(new_search_states);
            continue;
        }

        //  else: read is completely encapsulated within allele
        if (cache.empty) {
            cache.set(SearchState{
                    SA_Interval{sa_index, sa_index},
                    VariantSitePath{
                            VariantLocus{site_marker, allele_id}
                    },
                    VariantSitePath{},
                    SearchVariantSiteState::within_variant_site
            });
            continue;
        }

        VariantSitePath current_path = {VariantLocus{site_marker, allele_id}};
        bool cache_has_same_path = current_path == cache.search_state.traversed_path;
        if (cache_has_same_path) {
            cache.update_sa_interval_max(sa_index);
            continue;
        } else {
            cache.flush(new_search_states);
            cache.set(SearchState{
                    SA_Interval{sa_index, sa_index},
                    current_path,
                    VariantSitePath{},
                    SearchVariantSiteState::within_variant_site
            });
        }
    }
    cache.flush(new_search_states);
    return new_search_states;
}

SearchStates gram::handle_allele_encapsulated_states(const SearchStates &search_states,
                                                     const PRG_Info &prg_info) {
    SearchStates new_search_states = {};

    for (const auto &search_state: search_states) {
        bool has_path = not search_state.traversed_path.empty();
        if (has_path) {
            new_search_states.emplace_back(search_state);
            continue;
        }

        SearchStates split_search_states = handle_allele_encapsulated_state(search_state,
                                                                            prg_info);
        for (const auto &split_search_state: split_search_states)
            new_search_states.emplace_back(split_search_state);
    }
    return new_search_states;
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

SearchStates gram::search_read_backwards(const Pattern &read,
                                         const Pattern &kmer,
                                         const KmerIndex &kmer_index,
                                         const PRG_Info &prg_info) {
    // Test if kmer has been indexed
    bool kmer_in_index = kmer_index.find(kmer) != kmer_index.end();
    if (not kmer_in_index)
        return SearchStates{};

    // Test if kmer has been indexed, but has no search states in prg
    auto kmer_index_search_states = kmer_index.at(kmer);
    if (kmer_index_search_states.empty())
        return kmer_index_search_states;


    // Reverse iterator + skipping through indexed kmer in read
    auto read_begin = read.rbegin();
    std::advance(read_begin, kmer.size());

    SearchStates new_search_states = kmer_index_search_states;

    for (auto it = read_begin; it != read.rend(); ++it) { /// Iterates end to start of read
        const int_Base &pattern_char = *it;
        new_search_states = process_read_char_search_states(pattern_char,
                                                            new_search_states,
                                                            prg_info);
        // Test if no mapping found upon character extension
        auto read_not_mapped = new_search_states.empty();
        if (read_not_mapped)
            break;
    }

    if (!new_search_states.empty()) gram::set_allele_ids(new_search_states, prg_info);
    new_search_states = handle_allele_encapsulated_states(new_search_states, prg_info);
    return new_search_states;
}


/**
 * Backward search followed by check whether the extended searched pattern maps somewhere in the prg.
 */
SearchState search_fm_index_base_backwards(const int_Base &pattern_char,
                                           const uint64_t char_first_sa_index,
                                           const SearchState &search_state,
                                           const PRG_Info &prg_info) {
    auto next_sa_interval = base_next_sa_interval(pattern_char,
                                                  char_first_sa_index,
                                                  search_state.sa_interval,
                                                  prg_info);
    //  An 'invalid' SA interval (i,j) is defined by i-1=j, which occurs when the read no longer maps anywhere in the prg.
    auto valid_sa_interval = next_sa_interval.first - 1 != next_sa_interval.second;
    if (not valid_sa_interval) { // Create an empty, invalid search state.
        SearchState new_search_state;
        new_search_state.invalid = true;
        return new_search_state;
    }

    auto new_search_state = search_state;
    new_search_state.sa_interval.first = next_sa_interval.first;
    new_search_state.sa_interval.second = next_sa_interval.second;
    return new_search_state;
}

SearchStates gram::process_read_char_search_states(const int_Base &pattern_char,
                                                   const SearchStates &old_search_states,
                                                   const PRG_Info &prg_info) {
    //  Before extending backward search with next character, check for variant markers in the current SA intervals
    //  This is the v part of vBWT.
    auto post_markers_search_states = process_markers_search_states(old_search_states,
                                                                    prg_info);
    //  Regular backward searching
    auto new_search_states = search_base_backwards(pattern_char,
                                                   post_markers_search_states,
                                                   prg_info);
    return new_search_states;
}


SA_Interval gram::base_next_sa_interval(const Marker &next_char,
                                        const SA_Index &next_char_first_sa_index,
                                        const SA_Interval &current_sa_interval,
                                        const PRG_Info &prg_info) {
    const auto &current_sa_start = current_sa_interval.first;
    const auto &current_sa_end = current_sa_interval.second;

    SA_Index sa_start_offset;
    if (current_sa_start <= 0)
        sa_start_offset = 0;
    else {
        //  TODO: Consider deleting this if-clause, next_char should never be > 4, it probably never runs
        if (next_char > 4)
            sa_start_offset = prg_info.fm_index.bwt.rank(current_sa_start, next_char);
        else {
            sa_start_offset = dna_bwt_rank(current_sa_start,
                                           next_char,
                                           prg_info);
        }
    }

    SA_Index sa_end_offset;
    //  TODO: Consider deleting this if-clause, next_char should never be > 4, it probably never runs
    if (next_char > 4)
        sa_end_offset = prg_info.fm_index.bwt.rank(current_sa_end + 1, next_char);
    else {
        sa_end_offset = dna_bwt_rank(current_sa_end + 1,
                                     next_char,
                                     prg_info);
    }

    auto new_start = next_char_first_sa_index + sa_start_offset;
    auto new_end = next_char_first_sa_index + sa_end_offset - 1;
    return SA_Interval{new_start, new_end};
}

SearchStates gram::search_base_backwards(const int_Base &pattern_char,
                                         const SearchStates &search_states,
                                         const PRG_Info &prg_info) {
    // Compute the first occurrence of `pattern_char` in the suffix array. Necessary for backward search.
    auto char_alphabet_rank = prg_info.fm_index.char2comp[pattern_char];
    auto char_first_sa_index = prg_info.fm_index.C[char_alphabet_rank];

    SearchStates new_search_states = {};

    for (const auto &search_state: search_states) {
        SearchState new_search_state = search_fm_index_base_backwards(pattern_char,
                                                                      char_first_sa_index,
                                                                      search_state,
                                                                      prg_info);
        if (new_search_state.invalid)
            continue;
        new_search_states.emplace_back(new_search_state);
    }

    return new_search_states;
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



bool gram::marker_is_site_end(const Marker &allele_marker,
                              const SA_Index &sa_right_of_marker,
                              const PRG_Info &prg_info) {

    // Note: a possible alternative implementation is to use the fact that if the site ID of the
    // text position at `sa_right_of_marker` != `allele_marker` - 1,
    // then we are at the end of the site. (can use cov graph to find this site ID)
    auto marker_index = prg_info.fm_index[sa_right_of_marker] - 1;
    auto last_position = prg_info.last_allele_positions.at(allele_marker);
    bool is_site_end = last_position == marker_index;

    return is_site_end;
}

/**
 * Computes the full SA interval of a given allele marker.
 * Note that the way this is computed is robust to variant markers not being continuous:
 *      eg, could have site with markers 5/6 & another with 9/10 without a site with 7/8.
 *
 * `sigma`: the alphabet size
 * `char2comp`: maps a symbol of the alphabet to a positive integer in [0...sigma - 1]
 * `comp2char`: the inverse mapping of above
 */
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

AlleleId gram::get_allele_id(const SA_Index &allele_marker_sa_index,
                             const PRG_Info &prg_info) {
    //  What is the index of the character just before the marker allele in the original text?
    auto internal_allele_text_index = prg_info.fm_index[allele_marker_sa_index] - 1;
    auto allele_id = (AlleleId) prg_info.allele_mask[internal_allele_text_index];
    assert(allele_id > 0);
    return allele_id;
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


std::string gram::serialize_search_state(const SearchState &search_state) {
    std::stringstream ss;
    ss << "****** Search State ******" << std::endl;

    ss << "SA interval: ["
       << search_state.sa_interval.first
       << ", "
       << search_state.sa_interval.second
       << "]";
    ss << std::endl;

    if (not search_state.traversed_path.empty()) {
        ss << "Variant site path [marker, allele id]: " << std::endl;
        for (const auto &variant_site: search_state.traversed_path) {
            auto marker = variant_site.first;

            if (variant_site.second != 0) {
                const auto &allele_id = variant_site.second;
                ss << "[" << marker << ", " << allele_id << "]" << std::endl;
            }
        }
    }
    ss << "****** END Search State ******" << std::endl;
    return ss.str();
}


std::ostream &gram::operator<<(std::ostream &os, const SearchState &search_state) {
    os << serialize_search_state(search_state);
    return os;
}
