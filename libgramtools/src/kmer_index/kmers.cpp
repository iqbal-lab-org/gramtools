#include "kmer_index/kmers.hpp"


std::vector<PrgIndexRange> get_boundary_marker_indexes(const PRG_Info &prg_info) {
    std::vector<PrgIndexRange> boundary_marker_indexes;

    using MarkerIndex = uint64_t;
    std::unordered_map<Marker, MarkerIndex> start_indexes;

    for (uint64_t marker_count = 1;
         marker_count <= prg_info.markers_mask_count_set_bits;
         ++marker_count) {
        auto marker_index = prg_info.prg_markers_select(marker_count);
        uint64_t marker_char = prg_info.encoded_prg[marker_index];

        auto marker_is_site_boundary = marker_char % 2 != 0;
        if (not marker_is_site_boundary)
            continue;

        auto start_index_memorized = start_indexes.find(marker_char) != start_indexes.end();
        if (not start_index_memorized) {
            start_indexes[marker_char] = marker_index;
            continue;
        }

        auto start_index = start_indexes[marker_char];
        auto &end_index = marker_index;
        boundary_marker_indexes.emplace_back(PrgIndexRange {start_index, end_index});

        start_indexes.erase(marker_char);
    }
    return boundary_marker_indexes;
}


uint64_t find_site_end_boundary(const uint64_t &within_site_index,
                                const PRG_Info &prg_info) {
    auto last_prg_index = prg_info.encoded_prg.size() - 1;
    auto number_markers_before = prg_info.prg_markers_rank(within_site_index);

    for (uint64_t marker_count = number_markers_before + 1;
         marker_count <= prg_info.markers_mask_count_set_bits;
         ++marker_count) {
        auto marker_index = prg_info.prg_markers_select(marker_count);
        uint64_t marker_char = prg_info.encoded_prg[marker_index];

        auto marker_is_boundary = marker_char % 2 != 0;
        if (not marker_is_boundary)
            continue;

        auto char_is_last_in_prg = marker_index == last_prg_index;
        if (char_is_last_in_prg)
            return marker_index;

        auto next_char_within_allele = prg_info.allele_mask[marker_index + 1] != 0;
        if (next_char_within_allele)
            continue;
        return marker_index;
    }
}


uint64_t get_kmer_region_end_index(const uint64_t end_marker_index,
                                   const uint64_t max_read_size,
                                   const PRG_Info &prg_info) {
    auto last_prg_index = prg_info.encoded_prg.size() - 1;
    auto end_index = end_marker_index + max_read_size - 1;
    if (end_index > last_prg_index)
        end_index = last_prg_index;

    auto within_variant_site = prg_info.allele_mask[end_index] > 0
                               or prg_info.prg_markers_mask[end_index] != 0;
    if (within_variant_site) {
        end_index = find_site_end_boundary(end_index,
                                           prg_info);
    }
    return end_index;
}


std::vector<PrgIndexRange> get_kmer_region_ranges(std::vector<PrgIndexRange> &boundary_marker_indexes,
                                                  const uint64_t &max_read_size,
                                                  const PRG_Info &prg_info) {
    std::vector<PrgIndexRange> kmer_region_ranges;
    for (const auto &marker_indexes_range: boundary_marker_indexes) {
        auto &start_marker_index = marker_indexes_range.first;
        auto &end_marker_index = marker_indexes_range.second;

        auto kmer_region_start_index = start_marker_index;
        auto kmer_region_end_index = get_kmer_region_end_index(end_marker_index,
                                                               max_read_size,
                                                               prg_info);
        kmer_region_ranges.emplace_back(
                PrgIndexRange {kmer_region_start_index, kmer_region_end_index});
    }
    return kmer_region_ranges;
}


Patterns get_site_ordered_alleles(const uint64_t &within_site_index,
                                  const PRG_Info &prg_info) {
    auto site_end_index = find_site_end_boundary(within_site_index, prg_info);
    auto boundary_marker = prg_info.encoded_prg[site_end_index];

    int64_t current_index = (int64_t) site_end_index - 1;
    uint64_t current_char = 0;

    Patterns site_alleles;
    std::vector<Base> allele;

    while (current_char != boundary_marker and current_index >= 0) {
        current_char = prg_info.encoded_prg[current_index--];

        auto at_site_start_marker = current_char == boundary_marker;
        auto at_allele_marker = current_char > 4 and current_char % 2 == 0;
        if (at_site_start_marker or at_allele_marker) {
            std::reverse(allele.begin(), allele.end());
            site_alleles.push_back(allele);
            allele.clear();
            continue;
        }

        allele.emplace_back(current_char);
    }

    std::reverse(site_alleles.begin(), site_alleles.end());
    return site_alleles;
}


struct KmerTraversalState {
    uint64_t outside_site_start_index = 0;
    uint64_t last_marker_index = 0;
    uint64_t total_handled_sites_count = 0;
    uint64_t total_intersite_size = 0;
};


bool check_marker_in_kmer_range(const KmerTraversalState &traversal_state,
                                const uint64_t kmer_size) {
    auto kmer_distance_traversed = traversal_state.total_intersite_size
                                   + traversal_state.total_handled_sites_count;
    auto marker_in_range = kmer_distance_traversed + 1 <= kmer_size;
    return marker_in_range;
}


enum class MarkerHandlerStatus {
    marker_not_in_range,
    handled,
    unhandled,
};


MarkerHandlerStatus handle_allele_marker(KmerTraversalState &traversal_state,
                                         const uint64_t marker_index,
                                         const PRG_Info &prg_info) {
    auto marker_char = prg_info.encoded_prg[marker_index];

    auto at_allele_marker = marker_char % 2 == 0;
    if (not at_allele_marker)
        return MarkerHandlerStatus::unhandled;

    traversal_state.last_marker_index = marker_index;
    return MarkerHandlerStatus::handled;
}


MarkerHandlerStatus handle_first_marker_seen(std::list<uint64_t> &inrange_sites,
                                             KmerTraversalState &traversal_state,
                                             const uint64_t marker_index,
                                             const uint64_t kmer_size,
                                             const PRG_Info &prg_info) {
    auto at_first_marker = inrange_sites.empty();
    if (not at_first_marker)
        return MarkerHandlerStatus::unhandled;

    traversal_state.total_intersite_size = traversal_state.outside_site_start_index
                                           - marker_index;
    auto marker_in_range = check_marker_in_kmer_range(traversal_state,
                                                      kmer_size);

    if (not marker_in_range)
        return MarkerHandlerStatus::marker_not_in_range;

    traversal_state.last_marker_index = marker_index;
    inrange_sites.push_front(marker_index);
    return MarkerHandlerStatus::handled;
}


MarkerHandlerStatus handle_end_boundary_marker(std::list<uint64_t> &inrange_sites,
                                               KmerTraversalState &traversal_state,
                                               const uint64_t marker_index,
                                               const uint64_t kmer_size,
                                               const PRG_Info &prg_info) {
    auto marker_char = prg_info.encoded_prg[marker_index];
    auto at_boundary_marker = marker_char % 2 != 0;
    const auto &last_marker_index = traversal_state.last_marker_index;
    const auto &last_marker_char = prg_info.encoded_prg[last_marker_index];
    auto last_marker_was_boundary = last_marker_char % 2 != 0;
    auto at_boundary_end_marker = at_boundary_marker and last_marker_was_boundary;

    if (not at_boundary_end_marker)
        return MarkerHandlerStatus::unhandled;

    traversal_state.total_intersite_size += last_marker_index
                                            - marker_index
                                            - 1;
    auto marker_in_range = check_marker_in_kmer_range(traversal_state,
                                                      kmer_size);
    if (not marker_in_range)
        return MarkerHandlerStatus::marker_not_in_range;

    inrange_sites.push_front(marker_index);
    traversal_state.last_marker_index = marker_index;
    return MarkerHandlerStatus::handled;
}


MarkerHandlerStatus handle_start_boundary_marker(std::list<uint64_t> &inrange_sites,
                                                 KmerTraversalState &traversal_state,
                                                 const uint64_t marker_index,
                                                 const PRG_Info &prg_info) {
    auto marker_char = prg_info.encoded_prg[marker_index];
    auto at_boundary_marker = marker_char % 2 != 0;
    const auto &last_marker_index = traversal_state.last_marker_index;
    const auto &last_marker_char = prg_info.encoded_prg[last_marker_index];
    auto last_marker_was_boundary = last_marker_char % 2 != 0;
    auto at_boundary_start_marker = at_boundary_marker and not last_marker_was_boundary;

    if (not at_boundary_start_marker)
        return MarkerHandlerStatus::unhandled;

    traversal_state.total_handled_sites_count++;
    traversal_state.last_marker_index = marker_index;
    return MarkerHandlerStatus::handled;
}


MarkerHandlerStatus find_site_end_indexes(std::list<uint64_t> &inrange_sites,
                                          KmerTraversalState &traversal_state,
                                          const uint64_t marker_index,
                                          const uint64_t kmer_size,
                                          const PRG_Info &prg_info) {
    auto result_handled = [](const auto &result) {
        return result != MarkerHandlerStatus::unhandled;
    };

    auto result = handle_first_marker_seen(inrange_sites,
                                           traversal_state,
                                           marker_index,
                                           kmer_size,
                                           prg_info);
    if (result_handled(result))
        return result;

    result = handle_allele_marker(traversal_state,
                                  marker_index,
                                  prg_info);
    if (result_handled(result))
        return result;

    result = handle_end_boundary_marker(inrange_sites,
                                        traversal_state,
                                        marker_index,
                                        kmer_size,
                                        prg_info);
    if (result_handled(result))
        return result;

    result = handle_start_boundary_marker(inrange_sites, traversal_state, marker_index, prg_info);
    return result;
}


bool index_is_site_end_boundary(const uint64_t index,
                                const PRG_Info &prg_info) {
    auto at_last_prg_index = index == prg_info.encoded_prg.size() - 1;
    auto at_marker = prg_info.prg_markers_mask[index] == 1;
    bool at_site_end_boundary;
    if (not at_last_prg_index) {
        auto next_char_within_allele = prg_info.allele_mask[index + 1] > 0;
        at_site_end_boundary = at_marker and not next_char_within_allele;
    } else {
        at_site_end_boundary = at_marker;
    }
    return at_site_end_boundary;
}


std::list<uint64_t> sites_inrange_left(const uint64_t outside_site_start_index,
                                       const uint64_t kmer_size,
                                       const PRG_Info &prg_info) {
    const auto &start_index = outside_site_start_index;
    auto number_markers_before = prg_info.prg_markers_rank(start_index);
    auto at_site_end_boundary = index_is_site_end_boundary(start_index,
                                                           prg_info);
    if (at_site_end_boundary)
        number_markers_before += 1;

    std::list<uint64_t> inrange_sites;
    KmerTraversalState traversal_state;
    traversal_state.outside_site_start_index = outside_site_start_index;

    for (uint64_t marker_count = number_markers_before;
         marker_count > 0;
         --marker_count) {
        auto marker_index = prg_info.prg_markers_select(marker_count);
        auto result = find_site_end_indexes(inrange_sites,
                                            traversal_state,
                                            marker_index,
                                            kmer_size,
                                            prg_info);
        if (result == MarkerHandlerStatus::marker_not_in_range)
            break;
    }
    return inrange_sites;
}


std::pair<uint64_t, uint64_t> get_nonvariant_region(const uint64_t &site_end_boundary_index,
                                                    const PRG_Info &prg_info) {
    const auto &start_marker_index = site_end_boundary_index;
    auto last_prg_index = prg_info.encoded_prg.size() - 1;

    uint64_t nonvariant_region_start = 0;
    uint64_t nonvariant_region_end = 0;
    if (start_marker_index + 1 > last_prg_index)
        return std::make_pair(nonvariant_region_start, nonvariant_region_end);

    nonvariant_region_start = start_marker_index + 1;
    auto number_markers_before = prg_info.prg_markers_rank(start_marker_index);
    auto next_marker_offset = number_markers_before + 2;

    auto no_next_marker = next_marker_offset > prg_info.markers_mask_count_set_bits;
    if (no_next_marker)
        nonvariant_region_end = last_prg_index;
    else
        nonvariant_region_end = prg_info.prg_markers_select(next_marker_offset) - 1;
    return std::make_pair(nonvariant_region_start, nonvariant_region_end);
}


std::vector<Base> right_intersite_nonvariant_region(const uint64_t &site_end_boundary_index,
                                                    const PRG_Info &prg_info) {
    auto nonvariant_range = get_nonvariant_region(site_end_boundary_index,
                                                  prg_info);
    auto start = nonvariant_range.first;
    auto end = nonvariant_range.second;
    auto region_size = end - start + 1;

    std::vector<Base> nonvariant_region;
    nonvariant_region.reserve(region_size);

    for (auto i = start; i <= end; ++i) {
        auto base = (Base) prg_info.encoded_prg[i];
        nonvariant_region.push_back(base);
    }
    return nonvariant_region;
}


std::vector<Base> extract_simple_reverse_kmer(const uint64_t kmer_end_index,
                                              const uint64_t kmer_size,
                                              const PRG_Info &prg_info) {
    std::vector<Base> reverse_kmer;
    reverse_kmer.reserve(kmer_size);

    uint64_t kmer_start_index = 0;
    auto start_index_is_valid = (int64_t) kmer_end_index
                                - (int64_t) kmer_size + 1 >= 0;
    if (start_index_is_valid)
        kmer_start_index = kmer_end_index - kmer_size + 1;
    else
        return reverse_kmer;

    for (uint64_t i = kmer_end_index; i >= kmer_start_index; --i) {
        auto base = (Base) prg_info.encoded_prg[i];
        reverse_kmer.push_back(base);
        if (i == 0)
            break;
    }
    return reverse_kmer;
}


uint64_t find_site_start_boundary(const uint64_t &end_boundary_index,
                                  const PRG_Info &prg_info) {
    auto target_marker = prg_info.encoded_prg[end_boundary_index];

    uint64_t current_index = end_boundary_index;
    uint64_t current_marker = 0;
    auto number_markers_before = prg_info.prg_markers_rank(current_index);

    while (current_marker != target_marker) {
        current_index = prg_info.prg_markers_select(number_markers_before);
        current_marker = prg_info.encoded_prg[current_index];
        if (number_markers_before == 1)
            break;
        number_markers_before--;
    }
    return current_index;
}


Pattern get_pre_site_part(const uint64_t site_end_boundary,
                           const uint64_t kmer_size,
                           const PRG_Info &prg_info) {
    auto first_site_start_boundary = find_site_start_boundary(site_end_boundary,
                                                              prg_info);
    Pattern pre_site_part = {};
    if (first_site_start_boundary != 0) {
        int64_t end_index = first_site_start_boundary - kmer_size - 1;
        if (end_index < 0)
            end_index = 0;

        for (int64_t i = first_site_start_boundary - 1;
             i >= end_index; --i) {
            auto base = (Base) prg_info.encoded_prg[i];
            if (base > 4)
                break;
            pre_site_part.push_back(base);
        }
        std::reverse(pre_site_part.begin(), pre_site_part.end());
    }
    return pre_site_part;
}


void add_pre_site_region(std::list<Patterns> &region_parts,
                         const std::list<uint64_t> &inrange_sites,
                         const uint64_t kmer_size,
                         const PRG_Info &prg_info) {
    auto first_site_end_boundary = inrange_sites.front();
    Pattern pre_site_part = get_pre_site_part(first_site_end_boundary,
                                               kmer_size,
                                               prg_info);
    if (not pre_site_part.empty())
        region_parts.emplace_back(Patterns {pre_site_part});
}


void add_site_regions(std::list<Patterns> &region_parts,
                      const std::list<uint64_t> &inrange_sites,
                      const PRG_Info &prg_info) {
    auto site_count = 0;
    for (const auto &end_boundary_index: inrange_sites) {
        auto ordered_alleles = get_site_ordered_alleles(end_boundary_index,
                                                        prg_info);
        region_parts.push_back(ordered_alleles);

        auto at_last_site = site_count++ == inrange_sites.size() - 1;
        if (at_last_site)
            continue;

        auto nonvariant_region = right_intersite_nonvariant_region(end_boundary_index,
                                                                   prg_info);
        region_parts.emplace_back(Patterns {nonvariant_region});
    }
}


void add_post_site_regions(std::list<Patterns> &region_parts,
                           const std::list<uint64_t> &inrange_sites,
                           const uint64_t kmer_size,
                           const PRG_Info &prg_info) {
    auto end_boundary_index = inrange_sites.back();
    if (end_boundary_index == prg_info.encoded_prg.size() - 1)
        return;

    auto index = end_boundary_index + 1;
    uint64_t number_consumed_kmer_bases = 0;

    Pattern nonvariant_region = {};

    while (number_consumed_kmer_bases < kmer_size + 1
           and index <= prg_info.encoded_prg.size() - 1) {

        auto within_site = prg_info.allele_mask[index] > 0
                           or prg_info.prg_markers_mask[index] != 0;
        if (not within_site) {
            auto base = (Base) prg_info.encoded_prg[index];
            nonvariant_region.push_back(base);

            index++;
            number_consumed_kmer_bases++;
            continue;
        }

        if (not nonvariant_region.empty()) {
            region_parts.emplace_back(Patterns {nonvariant_region});
            nonvariant_region = {};
        }

        auto site_end_boundary = find_site_end_boundary(index, prg_info);
        auto ordered_alleles = get_site_ordered_alleles(site_end_boundary,
                                                        prg_info);
        region_parts.push_back(ordered_alleles);

        if (site_end_boundary == prg_info.encoded_prg.size() - 1)
            break;
        else
            index = site_end_boundary + 1;
        number_consumed_kmer_bases++;
    }

    if (not nonvariant_region.empty())
        region_parts.emplace_back(Patterns {nonvariant_region});
}


std::list<Patterns> get_kmer_size_region_parts(const uint64_t &current_range_end_index,
                                                    const std::list<uint64_t> &inrange_sites,
                                                    const uint64_t kmer_size,
                                                    const PRG_Info &prg_info) {
    std::list<Patterns> region_parts = {};
    add_pre_site_region(region_parts,
                        inrange_sites,
                        kmer_size,
                        prg_info);
    add_site_regions(region_parts,
                     inrange_sites,
                     prg_info);
    add_post_site_regions(region_parts,
                          inrange_sites,
                          kmer_size,
                          prg_info);
    return region_parts;
}


bool update_allele_index_path(std::vector<uint64_t> &allele_current_index,
                              const std::vector<uint64_t> &allele_counts) {
    bool update_possible_flag = false;

    // find rightmost index to increase
    int64_t i = (int64_t) allele_current_index.size() - 1;
    for (; i >= 0; --i) {
        auto at_last_valid_index = allele_current_index[i] + 1 == allele_counts[i];
        if (not at_last_valid_index) {
            update_possible_flag = true;
            break;
        }
    }

    if (not update_possible_flag)
        return update_possible_flag;

    allele_current_index[i]++;
    i++;

    for (; i < allele_current_index.size(); ++i)
        allele_current_index[i] = 0;
    return update_possible_flag;
}


Patterns get_paths_from_parts(const std::list<Patterns> &region_parts) {
    uint64_t number_of_paths_expected = 1;
    for (const auto &ordered_alleles: region_parts) {
        uint64_t number_of_alleles = ordered_alleles.size();
        number_of_paths_expected *= number_of_alleles;
    }

    std::vector<uint64_t> allele_current_index(region_parts.size(), 0);
    std::vector<uint64_t> allele_counts;
    allele_counts.reserve(region_parts.size());
    for (const auto &ordered_alleles: region_parts)
        allele_counts.push_back(ordered_alleles.size());

    Patterns paths;
    while (paths.size() < number_of_paths_expected) {
        Pattern path;

        uint64_t i = 0;
        for (const auto &ordered_alleles: region_parts) {
            auto allele_index = allele_current_index[i];
            const auto &allele = ordered_alleles[allele_index];
            path.insert(path.end(), allele.begin(), allele.end());
            i++;
        }

        paths.emplace_back(path);
        bool update_possible_flag = update_allele_index_path(allele_current_index,
                                                             allele_counts);
        if (not update_possible_flag)
            break;
    }
    return paths;
}


unordered_vector_set<Pattern> get_reverse_kmers_from_path(const Pattern &path,
                                                          const uint64_t &kmer_size) {
    unordered_vector_set<Pattern> reverse_kmers;
    for (int64_t i = path.size() - 1; i >= kmer_size - 1; --i) {
        Pattern reverse_kmer;
        reverse_kmer.reserve(kmer_size);
        for (int64_t j = i; j >= i - (int64_t) kmer_size + 1; --j)
            reverse_kmer.push_back(path[j]);
        reverse_kmers.insert(reverse_kmer);
    }
    return reverse_kmers;
}


unordered_vector_set<Pattern> extract_variant_reverse_kmers(uint64_t &current_range_end_index,
                                                             const std::list<uint64_t> &inrange_sites,
                                                             const uint64_t kmer_size,
                                                             const PRG_Info &prg_info) {
    auto region_parts = get_kmer_size_region_parts(current_range_end_index,
                                                   inrange_sites,
                                                   kmer_size,
                                                   prg_info);
    auto paths = get_paths_from_parts(region_parts);
    unordered_vector_set<Pattern> all_reverse_kmers;

    for (const auto &path: paths) {
        auto reverse_kmers = get_reverse_kmers_from_path(path, kmer_size);
        all_reverse_kmers.insert(reverse_kmers.begin(), reverse_kmers.end());
    }

    auto first_site_end_boundary = inrange_sites.front();
    auto first_site_start_boundary = find_site_start_boundary(first_site_end_boundary,
                                                              prg_info);
    if (first_site_start_boundary == 0)
        current_range_end_index = first_site_start_boundary;
    else
        current_range_end_index = first_site_start_boundary - 1;
    return all_reverse_kmers;
}


unordered_vector_set<Pattern> get_reverse_kmers_from_region(const PrgIndexRange &kmer_region_range,
                                                            const uint64_t &kmer_size,
                                                            const PRG_Info &prg_info) {
    const auto &region_start = kmer_region_range.first;
    const auto &region_end = kmer_region_range.second;

    unordered_vector_set<Pattern> all_reverse_kmers = {};

    for (auto current_index = region_end;
         current_index >= region_start;
         --current_index) {
        auto current_index_is_valid = (int64_t) current_index
                                      - (int64_t) kmer_size + 1 >= 0;
        if (not current_index_is_valid)
            break;

        auto inrange_sites = sites_inrange_left(current_index,
                                                kmer_size,
                                                prg_info);

        auto sites_in_range = not inrange_sites.empty();
        if (sites_in_range) {
            auto reverse_kmers = extract_variant_reverse_kmers(current_index,
                                                               inrange_sites,
                                                               kmer_size,
                                                               prg_info);
            all_reverse_kmers.insert(reverse_kmers.begin(), reverse_kmers.end());
            if (current_index == 0)
                break;
            continue;
        }

        auto within_site = prg_info.allele_mask[current_index] > 0
                           or prg_info.prg_markers_mask[current_index] != 0;
        if (not within_site) {
            auto reverse_kmer = extract_simple_reverse_kmer(current_index,
                                                            kmer_size,
                                                            prg_info);
            if (reverse_kmer.empty())
                break;
            all_reverse_kmers.insert(reverse_kmer);
            continue;
        }
    }
    return all_reverse_kmers;
}


std::vector<PrgIndexRange> combine_overlapping_regions(const std::vector<PrgIndexRange> &kmer_region_ranges) {
    auto sorted_ranges = kmer_region_ranges;

    auto ordering_condition = [](const auto &lhs, const auto &rhs) {
        if (lhs.first < rhs.first)
            return true;
        if (lhs.first == rhs.first)
            return lhs.second < rhs.second;
        return false;
    };

    std::sort(sorted_ranges.begin(),
              sorted_ranges.end(),
              ordering_condition);

    std::vector<PrgIndexRange> reduced_ranges;
    PrgIndexRange last_range_seen = {0, 0};
    for (const auto &range: sorted_ranges) {

        auto at_first_range = last_range_seen.first == 0
                              and last_range_seen.second == 0;
        if (at_first_range) {
            last_range_seen.first = range.first;
            last_range_seen.second = range.second;
            continue;
        }

        auto last_range_not_overlap = last_range_seen.second < range.first;
        if (last_range_not_overlap) {
            reduced_ranges.push_back(last_range_seen);
            last_range_seen.first = range.first;
            last_range_seen.second = range.second;
            continue;
        }

        auto range_completely_encapsulated = last_range_seen.second > range.second;
        if (range_completely_encapsulated)
            continue;

        last_range_seen.second = range.second;
    }

    auto no_last_range_recorded = last_range_seen.first == 0
                                  and last_range_seen.second == 0;
    if (no_last_range_recorded)
        return reduced_ranges;

    if (reduced_ranges.empty()) {
        reduced_ranges.push_back(last_range_seen);
        return reduced_ranges;
    }

    auto last_range_not_added = reduced_ranges.back() != last_range_seen;
    if (last_range_not_added)
        reduced_ranges.push_back(last_range_seen);
    return reduced_ranges;
}


ordered_vector_set<Pattern> get_all_reverse_kmers(const Parameters &parameters,
                                                   const PRG_Info &prg_info) {
    auto boundary_marker_indexes = get_boundary_marker_indexes(prg_info);
    auto kmer_region_ranges = get_kmer_region_ranges(boundary_marker_indexes,
                                                     parameters.max_read_size,
                                                     prg_info);
    kmer_region_ranges = combine_overlapping_regions(kmer_region_ranges);
    ordered_vector_set<Pattern> all_reverse_kmers;
    for (const auto &kmer_region_range: kmer_region_ranges) {
        auto reverse_kmers = get_reverse_kmers_from_region(kmer_region_range,
                                                           parameters.kmers_size,
                                                           prg_info);
        all_reverse_kmers.insert(reverse_kmers.begin(), reverse_kmers.end());
    }
    return all_reverse_kmers;
}


std::vector<Pattern> get_prefix_diffs(const std::vector<Pattern> &kmers) {
    std::vector<Pattern> prefix_diffs = {};
    Pattern last_full_kmer = {};
    for (const auto &kmer: kmers) {
        if (last_full_kmer.empty()) {
            last_full_kmer = kmer;
            prefix_diffs.push_back(last_full_kmer);
            continue;
        }

        bool prefix_found_flag = false;
        std::list<Base> prefix_diff_list = {};

        for (int64_t i = last_full_kmer.size() - 1; i >= 0; --i) {
            auto &base = kmer[i];
            auto &last_full_base = last_full_kmer[i];

            if (base != last_full_base)
                prefix_found_flag = true;

            if (prefix_found_flag)
                prefix_diff_list.push_front(base);
        }
        last_full_kmer = kmer;
        Pattern prefix_diff{std::make_move_iterator(std::begin(prefix_diff_list)),
                            std::make_move_iterator(std::end(prefix_diff_list))};
        prefix_diffs.push_back(prefix_diff);
    }
    return prefix_diffs;
}


std::vector<Pattern> reverse_kmers_inplace(const ordered_vector_set<Pattern> &reverse_kmers) {
    std::vector<Pattern> kmers;
    for (auto reverse_kmer: reverse_kmers) {
        std::reverse(reverse_kmer.begin(), reverse_kmer.end());
        auto &kmer = reverse_kmer;
        kmers.emplace_back(kmer);
    }
    return kmers;
}


std::vector<Pattern> get_all_ordered_kmers(const Parameters &parameters,
                                            const PRG_Info &prg_info) {
    auto reverse_kmers = get_all_reverse_kmers(parameters,
                                               prg_info);
    auto kmers = reverse_kmers_inplace(reverse_kmers);
    return kmers;
}


std::vector<Pattern> get_kmer_prefix_diffs(const Parameters &parameters,
                                            const PRG_Info &prg_info) {
    auto kmers = get_all_ordered_kmers(parameters,
                                       prg_info);
    auto prefix_diffs = get_prefix_diffs(kmers);
    return prefix_diffs;
}
