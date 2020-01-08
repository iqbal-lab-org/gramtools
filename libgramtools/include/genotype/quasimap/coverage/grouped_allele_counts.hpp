/** @file
* Defines coverage related operations for base-level allele coverage.
*/
#include "genotype/quasimap/search_types.hpp"
#include "genotype/quasimap/coverage/types.hpp"
#include "common.hpp"


#ifndef GRAMTOOLS_GROUPED_ALLELE_COUNTS_HPP
#define GRAMTOOLS_GROUPED_ALLELE_COUNTS_HPP

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
                                       const Parameters &parameters);
        }
    }

    /**
     * Assigns a unique group ID to each distinct `gram::AlleleIds` group.
     */
    AlleleGroupHash hash_allele_groups(const SitesGroupedAlleleCounts &sites);

    /**
     * String-serialise a single site count.
     * Outputs an allele group ID and a count of reads mapped to that allele ID combination.
     * If no read has mapped to the site, outputs an empty entry ("{}").
     */
    std::string dump_site(const AlleleGroupHash &allele_ids_groups_hash,
                          const GroupedAlleleCounts &site);

    /**
     * String-serialise site counts in JSON format.
     * Site counts is an array where each element refers to a site.
     * @see dump_site()
     */
    std::string dump_site_counts(const AlleleGroupHash &allele_ids_groups_hash,
                                 const SitesGroupedAlleleCounts &sites);

    std::string dump_allele_groups(const AlleleGroupHash &allele_ids_groups_hash);

    std::string dump_grouped_allele_counts(const SitesGroupedAlleleCounts &sites);
}

#endif //GRAMTOOLS_GROUPED_ALLELE_COUNTS_HPP
