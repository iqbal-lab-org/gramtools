/**@file
 * Defines coverage related operations for base-level allele coverage.
 */
#include "quasimap/search_types.hpp"
#include "quasimap/coverage/types.hpp"
#include <optional>


#ifndef GRAMTOOLS_ALLELE_BASE_HPP
#define GRAMTOOLS_ALLELE_BASE_HPP

namespace gram {

    namespace coverage {
        namespace generate {
            /**
             * Produce base-level coverage recording structure.
             * @see types.hpp
             */
            SitesAlleleBaseCoverage allele_base_structure(const PRG_Info &prg_info);
        }

        namespace per_base{

            class Traverser{
            public:
                Traverser() {;}
                Traverser(node_access start_point, VariantSitePath traversed_loci, std::size_t read_size);

                std::optional<covG_ptr> next_Node();
                /*
                 * Getters
                 */
                std::pair<uint32_t, uint32_t> get_node_interval(){
                    return {start_pos, end_pos};
                }
                std::size_t get_remaining_bases(){
                    return bases_remaining;
                }

                void go_to_next_site();
                void move_past_single_edge_node();
                void update_coordinates();
                void assign_end_position();
                void choose_allele();

            private:
                covG_ptr cur_Node;
                std::size_t bases_remaining;
                VariantSitePath traversed_loci;
                uint32_t traversed_index;
                bool first_node;
                uint32_t start_pos;
                uint32_t end_pos;
            };
        }

        namespace record {
            /**
             * Record base-level coverage for selected `SearchStates`.
             * `SearchStates`, can have different mapping instances going through the same `VariantLocus`.
             * The `SitesCoverageBoundaries` structure avoids recording the same base more than once in that case.
             */
            void allele_base(Coverage &coverage,
                             const SearchStates &search_states,
                             const uint64_t &read_length,
                             const PRG_Info &prg_info);
        }

        namespace dump {
            /**
             * String serialise the coverage information in JSON format and write it to disk.
             */
            void allele_base(const Coverage &coverage,
                             const Parameters &parameters);
        }
    }

    std::string dump_allele_base_coverage(const SitesAlleleBaseCoverage &sites);

    /**
     * Compute the (start,end) positions in the prg of a variant site marker.
     * @param site_marker the variant site marker to find the positions of.
     */
    std::pair<uint64_t, uint64_t> site_marker_prg_indexes(const uint64_t &site_marker, const PRG_Info &prg_info);

    using SitesCoverageBoundaries = PairHashMap<VariantLocus, uint64_t>; /**< For a given `VariantLocus`, gives the last allele base position recorded.*/

    /**
     * Increments each traversed base's coverage in traversed allele.
     * @return The number of bases of the read processed forwards.
     */
    uint64_t set_site_base_coverage(Coverage &coverage,
                                    SitesCoverageBoundaries &sites_coverage_boundaries,
                                    const VariantLocus &path_element,
                                    const uint64_t allele_coverage_offset,
                                    const uint64_t max_bases_to_set);

    /**
    * Computes the difference between an index into an allele and the index of the allele's start.
    */
    uint64_t allele_start_offset_index(const uint64_t within_allele_prg_index, const PRG_Info &prg_info);

}

#endif //GRAMTOOLS_ALLELE_BASE_HPP
