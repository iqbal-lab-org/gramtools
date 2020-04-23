/** @file
* Defines coverage related operations for base-level allele coverage.
*/


#ifndef GRAMTOOLS_GROUPED_ALLELE_COUNTS_HPP
#define GRAMTOOLS_GROUPED_ALLELE_COUNTS_HPP

#include <nlohmann/json.hpp>

#include "genotype/quasimap/search/types.hpp"
#include "genotype/quasimap/coverage/types.hpp"
#include "coverage_common.hpp"

using JSON = nlohmann::json;

namespace gram {
    namespace coverage {
        namespace generate {
            /** Sets up the structure for recording grouped allele counts.
             * The structure is a vector holding, for each variant site of the prg, an unordered_map that can
             * associate together alleles mapped by the same read.
             * @see SitesGroupedAlleleCounts
             */
            SitesGroupedAlleleCounts grouped_allele_counts(const PRG_Info &prg_info);
        }

        namespace record {
            /**
             * Records allele group counts per site.
             * @see GroupedAlleleCounts
             * @note Single alleles also get registered as 'groups'.
             */
            void grouped_allele_counts(Coverage &coverage,
                                       uniqueLoci const& compatible_loci);
        }

        namespace dump {
            /**
             * Write grouped allele coverage to disk in JSON format.
             */
            void grouped_allele_counts(const Coverage &coverage,
                                       const GenotypeParams &parameters);
        }
    }

    /**
     * Assigns a unique group ID to each distinct `gram::AlleleIds` group.
     */
    AlleleGroupHash hash_allele_groups(const SitesGroupedAlleleCounts &sites);

    struct numericStringComparator {
        bool operator()(const std::string& lhs, const std::string& rhs) const {
            return std::stoi(lhs) < std::stoi(rhs);
        }
    };

    using GroupIDToCounts = std::map<std::string, CovCount, numericStringComparator>;
    using SitesGroupIDToCounts = std::vector<GroupIDToCounts>;
    using GroupIDToAlleles = std::map<std::string, AlleleIds, numericStringComparator>;

    SitesGroupIDToCounts get_group_id_counts(SitesGroupedAlleleCounts const& sites,
                                             AlleleGroupHash const& allele_ids_groups_hash);

    GroupIDToAlleles get_group_id_alleles(AlleleGroupHash const& allele_ids_group_hash);


    JSON get_json(const SitesGroupedAlleleCounts &sites, AlleleGroupHash const &allele_ids_groups_hash);
}

#endif //GRAMTOOLS_GROUPED_ALLELE_COUNTS_HPP
