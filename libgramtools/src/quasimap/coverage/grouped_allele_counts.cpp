#include <fstream>
#include <vector>

#include "search.hpp"

#include "quasimap/utils.hpp"
#include "quasimap/coverage/grouped_allele_counts.hpp"


SitesGroupedAlleleCounts coverage::generate::grouped_allele_counts(const PRG_Info &prg_info) {
    uint64_t numer_of_variant_sites = get_number_of_variant_sites(prg_info);
    SitesGroupedAlleleCounts grouped_allele_counts(numer_of_variant_sites);
    return grouped_allele_counts;
}


void coverage::record::grouped_allele_counts(Coverage &coverage,
                                             const SearchStates &search_states) {
    using AlleleIdSet = std::set<AlleleId>;
    std::unordered_map<Marker, AlleleIdSet> site_allele_group;

    for (const auto &search_state: search_states) {
        for (const auto &variant_site: search_state.variant_site_path) {
            auto site_marker = variant_site.first;
            auto allell_id = variant_site.second;
            site_allele_group[site_marker].insert(allell_id);
        }
    }

    for (const auto &entry: site_allele_group) {
        auto site_marker = entry.first;
        auto allele_ids_set = entry.second;
        AlleleIds allele_ids;
        std::copy(allele_ids_set.begin(), allele_ids_set.end(),
                  std::back_inserter(allele_ids));

        auto min_boundary_marker = 5;
        auto site_coverage_index = (site_marker - min_boundary_marker) / 2;

        auto &site_coverage = coverage.grouped_allele_counts[site_coverage_index];
        site_coverage[allele_ids] += 1;
    }
}


AlleleGroupHash hash_allele_groups(const SitesGroupedAlleleCounts &sites) {
    AlleleGroupHash allele_ids_groups_hash;
    uint64_t group_hash = 0;
    for (const auto &site: sites) {
        for (const auto &allele_group: site) {
            auto allele_ids_group = allele_group.first;
            auto group_seen = allele_ids_groups_hash.find(allele_ids_group)
                              != allele_ids_groups_hash.end();
            if (group_seen)
                continue;
            allele_ids_groups_hash[allele_ids_group] = group_hash++;
        }
    }
    return allele_ids_groups_hash;
}


std::string dump_site(const AlleleGroupHash &allele_ids_groups_hash,
                      const GroupedAlleleCounts &site) {
    std::stringstream stream;
    stream << "{";
    auto i = 0;
    for (const auto &allele_entry: site) {
        auto allele_ids_group = allele_entry.first;
        uint64_t group_hash = allele_ids_groups_hash.at(allele_ids_group);
        auto count = allele_entry.second;

        stream << "\"" << (int) group_hash << "\":" << (int) count;
        if (i++ < site.size() - 1)
            stream << ",";
    }
    stream << "}";
    return stream.str();
}


std::string dump_site_counts(const AlleleGroupHash &allele_ids_groups_hash,
                             const SitesGroupedAlleleCounts &sites) {
    std::stringstream stream;
    stream << "\"site_counts\":[";
    auto i = 0;
    for (const auto &site: sites) {
        stream << dump_site(allele_ids_groups_hash, site);
        if (i++ < sites.size() - 1)
            stream << ",";
    }
    stream << "]";
    return stream.str();
}


std::string dump_allele_groups(const AlleleGroupHash &allele_ids_groups_hash) {
    std::stringstream stream;
    stream << "\"allele_groups\":{";
    auto i = 0;
    for (const auto &entry: allele_ids_groups_hash) {
        auto allele_ids_group = entry.first;
        uint64_t group_hash = entry.second;
        stream << "\"" << (int) group_hash << "\":[";
        auto j = 0;
        for (const auto &allele_id: allele_ids_group) {
            stream << (int) allele_id;
            if (j++ < allele_ids_group.size() - 1)
                stream << ",";
        }
        stream << "]";
        if (i++ < allele_ids_groups_hash.size() - 1)
            stream << ",";
    }
    stream << "}";
    return stream.str();
}


std::string dump_grouped_allele_counts(const SitesGroupedAlleleCounts &sites) {
    auto allele_ids_groups_hash = hash_allele_groups(sites);
    std::stringstream stream;
    stream << "{\"grouped_allele_counts\":{";
    stream << dump_site_counts(allele_ids_groups_hash, sites);
    stream << ",";
    stream << dump_allele_groups(allele_ids_groups_hash);
    stream << "}}";
    return stream.str();
}


void coverage::dump::grouped_allele_counts(const Coverage &coverage,
                                           const Parameters &parameters) {
    std::string json_string = dump_grouped_allele_counts(coverage.grouped_allele_counts);
    std::ofstream file;
    file.open(parameters.grouped_allele_counts_fpath);
    file << json_string << std::endl;
}