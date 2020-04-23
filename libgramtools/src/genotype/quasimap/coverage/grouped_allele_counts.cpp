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
    std::unordered_map<Marker, AlleleIdSet> site_allele_group;

    // Loop through all loci and record the alleles compatible with each site
    for (const auto &locus : compatible_loci) {
            auto site_marker = locus.first;
            auto allele_id = locus.second;
            site_allele_group[site_marker].insert(allele_id);
        }

    // Loop through the variant site markers traversed at least once by the read.
    for (const auto &entry: site_allele_group) {
        auto site_marker = entry.first;
        auto allele_ids_set = entry.second;
        AlleleIds allele_ids;
        std::copy(allele_ids_set.begin(), allele_ids_set.end(),
                  std::back_inserter(allele_ids));

        auto site_index = siteID_to_index(site_marker);

        // Get the map between allele Ids and counts.
        auto &site_coverage = coverage.grouped_allele_counts[site_index];
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

SitesGroupIDToCounts gram::get_group_id_counts(SitesGroupedAlleleCounts const& sites,
                                               AlleleGroupHash const& allele_ids_groups_hash){
    SitesGroupIDToCounts result(sites.size());
    std::size_t num_processed{0};
    for (auto const& site: sites){
        GroupIDToCounts site_groups;
        for (auto const& equiv_class_count: site){
           auto group_id = allele_ids_groups_hash.at(equiv_class_count.first);
           site_groups[std::to_string(group_id)] = equiv_class_count.second;
        }
        result.at(num_processed++) = site_groups;
    }
    return result;
}

GroupIDToAlleles gram::get_group_id_alleles(AlleleGroupHash const& allele_ids_group_hash){
    GroupIDToAlleles result;
    for (auto const& entry: allele_ids_group_hash)
        result[std::to_string(entry.second)] = entry.first;
    return result;
}

JSON gram::get_json(const SitesGroupedAlleleCounts &sites, AlleleGroupHash const &allele_ids_groups_hash) {
    auto group_id_alleles = get_group_id_alleles(allele_ids_groups_hash);
    auto group_id_counts = get_group_id_counts(sites, allele_ids_groups_hash);
    JSON result;
    result["grouped_allele_counts"]["site_counts"] = group_id_counts;
    result["grouped_allele_counts"]["allele_groups"] = group_id_alleles;
    return result;
}


void coverage::dump::grouped_allele_counts(const Coverage &coverage,
                                           const GenotypeParams &parameters) {
    auto allele_ids_groups_hash = hash_allele_groups(coverage.grouped_allele_counts);
    auto json = get_json(coverage.grouped_allele_counts, allele_ids_groups_hash);
    std::ofstream file;
    file.open(parameters.grouped_allele_counts_fpath);
    file << json << std::endl;
}