#include <fstream>
#include <vector>

#include "genotype/quasimap/coverage/grouped_allele_counts.hpp"


using namespace gram;


SitesGroupedAlleleCounts coverage::generate::grouped_allele_counts(const PRG_Info &prg_info) {
    uint64_t number_of_variant_sites = prg_info.num_variant_sites;
    // Stores as many empty maps (associating allele IDs to mapped read counts) as there are variant sites in the prg.
    SitesGroupedAlleleCounts grouped_allele_counts(number_of_variant_sites);
    return grouped_allele_counts;
}


void coverage::record::grouped_allele_counts(Coverage &coverage,
                                             uniqueLoci const& compatible_loci) {
    // We will store, for each variant site `Marker`, which alleles are traversed across
    // **all** (selected, ie site-equivalent) mapping instances of the processed read.
    using AlleleIdSet = std::set<AlleleId>;
    std::unordered_map<Marker, AlleleIdSet> site_allele_group;

    // Loop through all loci and record the alleles compatible with each site
    for (const auto &locus : compatible_loci) {
            auto site_marker = locus.first;
            auto allele_id = locus.second - 1;
            site_allele_group[site_marker].insert(allele_id);
        }

    // Loop through the variant site markers traversed at least once by the read.
    for (const auto &entry: site_allele_group) {
        auto site_marker = entry.first;
        auto allele_ids_set = entry.second;
        AlleleIds allele_ids;
        std::copy(allele_ids_set.begin(), allele_ids_set.end(),
                  std::back_inserter(allele_ids));

        auto min_boundary_marker = 5;
        // Which site entry in `grouped_allele_counts` is concerned.
        auto site_coverage_index = (site_marker - min_boundary_marker) / 2;

        // Get the map between allele Ids and counts.
        auto &site_coverage = coverage.grouped_allele_counts[site_coverage_index];
        #pragma omp critical
        // Note: if the key does not already exists, creates a key value pair **and** initialises the value to 0.
        site_coverage[allele_ids] += 1;
    }
}


AlleleGroupHash gram::hash_allele_groups(const SitesGroupedAlleleCounts &sites) {
    AlleleGroupHash allele_ids_groups_hash;
    uint64_t group_ID = 0;
    // Loop through all allele id groups across all variant sites.
    for (const auto &site: sites) {
        for (const auto &allele_group: site) {
            auto allele_ids_group = allele_group.first;
            auto group_seen = allele_ids_groups_hash.find(allele_ids_group)
                              != allele_ids_groups_hash.end();
            // Group already has an ID, continue.
            if (group_seen)
                continue;
            allele_ids_groups_hash[allele_ids_group] = group_ID++;
        }
    }
    return allele_ids_groups_hash;
}


std::string gram::dump_site(const AlleleGroupHash &allele_ids_groups_hash,
                            const GroupedAlleleCounts &site) {
    std::stringstream stream;
    stream << "{";
    auto i = 0;
    for (const auto &allele_entry: site) {
        auto allele_ids_group = allele_entry.first;
        uint64_t group_ID = allele_ids_groups_hash.at(allele_ids_group);
        auto count = allele_entry.second;

        stream << "\"" << (int) group_ID << "\":" << (int) count;
        if (i++ < site.size() - 1)
            stream << ",";
    }
    stream << "}";
    return stream.str();
}


std::string gram::dump_site_counts(const AlleleGroupHash &allele_ids_groups_hash,
                                   const SitesGroupedAlleleCounts &sites) {
    std::stringstream stream;
    stream << "\"site_counts\":[";
    auto i = 0;
    // Call dump_site() for each site, regardless of whether it has any coverage information at all.
    for (const auto &site: sites) {
        stream << dump_site(allele_ids_groups_hash, site);
        if (i++ < sites.size() - 1)
            stream << ",";
    }
    stream << "]";
    return stream.str();
}


std::string gram::dump_allele_groups(const AlleleGroupHash &allele_ids_groups_hash) {
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


std::string gram::dump_grouped_allele_counts(const SitesGroupedAlleleCounts &sites) {
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
                                           const GenotypeParams &parameters) {
    std::string json_string = dump_grouped_allele_counts(coverage.grouped_allele_counts);
    std::ofstream file;
    file.open(parameters.grouped_allele_counts_fpath);
    file << json_string << std::endl;
}