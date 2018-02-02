#include <cassert>
#include <fstream>
#include <vector>

#include "search/search.hpp"

#include "quasimap/utils.hpp"
#include "quasimap/coverage/allele_base.hpp"


SitesAlleleBaseCoverage coverage::generate::allele_base_structure(const PRG_Info &prg_info) {
    uint64_t numer_of_variant_sites = get_number_of_variant_sites(prg_info);
    SitesAlleleBaseCoverage allele_base_coverage(numer_of_variant_sites);

    const auto min_boundary_marker = 5;

    uint64_t allele_size = 0;
    Marker last_marker = 0;

    for (const auto &mask_value: prg_info.sites_mask) {
        auto within_allele = mask_value != 0;
        if (within_allele) {
            allele_size += 1;
            last_marker = mask_value;
            continue;
        }

        auto no_allele_to_flush = allele_size == 0;
        if (no_allele_to_flush)
            continue;

        BaseCoverage bases(allele_size);
        uint64_t variant_site_cover_index = (last_marker - min_boundary_marker) / 2;
        allele_base_coverage.at(variant_site_cover_index).emplace_back(bases);
        allele_size = 0;
    }
    return allele_base_coverage;
}


uint64_t allele_start_offset_index(const uint64_t within_allele_prg_index, const PRG_Info &prg_info) {
    // find the nearest left marker with rank and select
    auto number_markers_before = prg_info.prg_markers_rank(within_allele_prg_index);
    auto marker_index = prg_info.prg_markers_select(number_markers_before);
    auto offset = within_allele_prg_index - marker_index - 1;
    return offset;
}


uint64_t set_site_base_coverage(Coverage &coverage,
                                SitesCoverageBoundaries &sites_coverage_boundaries,
                                const VariantSite &path_element,
                                const uint64_t allele_coverage_offset,
                                const uint64_t max_set_bases) {
    auto marker = path_element.first;
    auto min_boundary_marker = 5;
    auto variant_site_coverage_index = (marker - min_boundary_marker) / 2;
    auto &site_coverage = coverage.allele_base_coverage.at(variant_site_coverage_index);

    auto allell_id = path_element.second;
    auto allele_coverage_index = allell_id - 1;
    auto &allele_coverage = site_coverage.at(allele_coverage_index);

    uint64_t index_end_boundary = std::min(max_set_bases, allele_coverage.size());
    uint64_t count_bases_consumed = index_end_boundary - allele_coverage_offset;

    uint64_t index_start_boundary = allele_coverage_offset;
    bool site_seen_previously = sites_coverage_boundaries.find(path_element) != sites_coverage_boundaries.end();
    if (site_seen_previously)
        index_start_boundary = std::max(allele_coverage_offset, sites_coverage_boundaries[path_element]);
    sites_coverage_boundaries[path_element] = index_end_boundary;

    for (uint64_t i = index_start_boundary; i < index_end_boundary; ++i)
        ++allele_coverage[i];
    return count_bases_consumed;
}


std::pair<uint64_t, uint64_t> site_marker_prg_indexes(const uint64_t &site_marker, const PRG_Info &prg_info) {
    auto alphabet_rank = prg_info.fm_index.char2comp[site_marker];
    auto first_sa_index = prg_info.fm_index.C[alphabet_rank];
    auto second_sa_index = first_sa_index + 1;

    auto first_prg_index = prg_info.fm_index[first_sa_index];
    auto second_prg_index = prg_info.fm_index[second_sa_index];

    if (first_prg_index < second_prg_index)
        return std::make_pair(first_prg_index, second_prg_index);
    else
        return std::make_pair(second_prg_index, first_prg_index);
}


uint64_t inter_site_base_count(const uint64_t &first_site_marker,
                               const uint64_t &second_site_marker,
                               const PRG_Info &prg_info) {
    auto first_site_prg_start_end = site_marker_prg_indexes(first_site_marker, prg_info);
    auto first_site_prg_end = first_site_prg_start_end.second;

    auto second_site_prg_start_end = site_marker_prg_indexes(second_site_marker, prg_info);
    auto second_site_prg_start = second_site_prg_start_end.first;

    return second_site_prg_start - first_site_prg_end - 1;
}


void sa_index_allele_base_coverage(Coverage &coverage,
                                   SitesCoverageBoundaries &sites_coverage_boundaries,
                                   const uint64_t &sa_index,
                                   const uint64_t &read_length,
                                   const SearchState &search_state,
                                   const PRG_Info &prg_info) {
    uint64_t read_bases_consumed = 0;
    uint64_t last_site_marker = 0;
    auto path_it = search_state.variant_site_path.begin();

    auto read_start_index = prg_info.fm_index[sa_index];
    auto start_site_marker = prg_info.sites_mask[read_start_index];
    bool read_starts_within_site = start_site_marker != 0;
    if (read_starts_within_site) {
        const auto &path_element = *path_it;
        last_site_marker = path_element.first;

        auto allele_coverage_offset = allele_start_offset_index(read_start_index, prg_info);
        auto max_set_bases = read_length - read_bases_consumed;
        read_bases_consumed += set_site_base_coverage(coverage,
                                                      sites_coverage_boundaries,
                                                      path_element,
                                                      allele_coverage_offset,
                                                      max_set_bases);
        ++path_it;
    } else {
        // forward to first path element
        const auto &path_element = *path_it;
        auto site_marker = path_element.first;
        auto site_prg_start_end = site_marker_prg_indexes(site_marker, prg_info);
        read_bases_consumed += site_prg_start_end.first - read_start_index;
    }

    auto last_path_it = search_state.variant_site_path.end();
    while (read_bases_consumed < read_length and path_it != last_path_it) {
        // take account of bases between sites
        const auto &path_element = *path_it;
        auto site_marker = path_element.first;

        if (last_site_marker != 0)
            read_bases_consumed += inter_site_base_count(last_site_marker, site_marker, prg_info);
        last_site_marker = site_marker;

        uint64_t allele_coverage_offset = 0;
        auto max_set_bases = read_length - read_bases_consumed;
        read_bases_consumed += set_site_base_coverage(coverage,
                                                      sites_coverage_boundaries,
                                                      path_element,
                                                      allele_coverage_offset,
                                                      max_set_bases);
        ++path_it;
    }
}


void coverage::record::allele_base(Coverage &coverage,
                                   const SearchStates &search_states,
                                   const uint64_t &read_length,
                                   const PRG_Info &prg_info) {
    SitesCoverageBoundaries sites_coverage_boundaries;

    for (const auto &search_state: search_states) {
        if (search_state.variant_site_path.empty())
            continue;

        auto first_sa_index = search_state.sa_interval.first;
        auto last_sa_index = search_state.sa_interval.second;
        for (auto sa_index = first_sa_index; sa_index <= last_sa_index; ++sa_index)
            sa_index_allele_base_coverage(coverage,
                                          sites_coverage_boundaries,
                                          sa_index,
                                          read_length,
                                          search_state,
                                          prg_info);
    }
}


std::string dump_allele(const BaseCoverage &allele) {
    std::stringstream stream;
    stream << "[";
    auto i = 0;
    for (const auto &base_coverage: allele) {
        stream << (int) base_coverage;
        if (i++ < allele.size() - 1)
            stream << ",";
    }
    stream << "]";
    return stream.str();
}


std::string dump_site(const AlleleCoverage &site) {
    std::stringstream stream;
    auto i = 0;
    for (const auto &allele: site) {
        stream << dump_allele(allele);
        if (i++ < site.size() - 1)
            stream << ",";
    }
    return stream.str();
}


std::string dump_sites(const SitesAlleleBaseCoverage &sites) {
    std::stringstream stream;
    auto i = 0;
    for (const auto &site: sites) {
        stream << "[";
        stream << dump_site(site);
        stream << "]";
        if (i++ < sites.size() - 1)
            stream << ",";
    }
    return stream.str();
}


std::string dump_allele_base_coverage(const SitesAlleleBaseCoverage &sites) {
    std::stringstream stream;
    stream << "{\"allele_base_counts\":[";
    stream << dump_sites(sites);
    stream << "]}";
    return stream.str();
}


void coverage::dump::allele_base(const Coverage &coverage,
                                 const Parameters &parameters) {
    std::string json_string = dump_allele_base_coverage(coverage.allele_base_coverage);
    std::ofstream file;
    file.open(parameters.allele_base_coverage_fpath);
    file << json_string << std::endl;
}