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

namespace gram {
    namespace coverage {
        namespace per_base {
            using node_coordinate = uint32_t;
            using node_coordinates = std::pair<node_coordinate, node_coordinate>;

        class InconsistentCovNodeCoordinates : public std::exception {
            public:
                InconsistentCovNodeCoordinates(std::string msg) : msg(msg) {;}
                const char * what() const throw(){
                    return msg.c_str();
                }
            private:
                std::string msg;
            };

            class DummyCovNode{
            public:
                DummyCovNode() = default;
                DummyCovNode(node_coordinate start_pos, node_coordinate end_pos, std::size_t node_size);
                bool operator==(DummyCovNode const& other) const{
                    return (node_size == other.node_size &&
                            start_pos == other.start_pos &&
                            end_pos == other.end_pos &&
                            node_size == other.node_size);
                }

                void extend_coordinates(node_coordinates coords);
                node_coordinates get_coordinates() const{return node_coordinates{start_pos, end_pos};}
            private:
                bool full;
                node_coordinate start_pos;
                node_coordinate end_pos;
                std::size_t node_size;
            };

            using realCov_to_dummyCov = std::map<covG_ptr, DummyCovNode>;

            /**
             * Class which produces all coverage node from the coverage graph that are in variant sites.
             * The choice of nodes at fork points is made using the set of `VariantLocus` traversed by a `SearchState`.
             *
             * Note the current assumption must be true: each node in a bubble has outdegree 1.
             * This is enforced in the `coverage_Graph` by having site boundary nodes flanking each bubble.
             */
            class Traverser {
            public:
                Traverser() {}

                Traverser(node_access start_point, VariantSitePath traversed_loci, std::size_t read_size);

                std::optional<covG_ptr> next_Node();

                /*
                 * Getters
                 */
                node_coordinates get_node_coordinates() {
                    return {start_pos, end_pos};
                }

                std::size_t get_remaining_bases() {
                    return bases_remaining;
                }

                /**
                 * Advances past all nodes with outdegree one, until we either:
                 *  - Find a node with outdegree > 1, so we choose an allelic node
                 *  - Consume all bases, so the traversal has ended.
                 */
                void go_to_next_site();

                /**
                 * First node gets special treatment.
                 * We can either start:
                 *  - Outside of a bubble: in which case, we move to the next node in a bubble
                 *  - In a bubble: in which case, we only call update_coordinates()
                 */
                void process_first_node();

                /**
                 * Consumes bases in the current node, and sets start and end coordinates.
                 * The start and end coordinates signal how much coverage should be recorded.
                 */
                void update_coordinates();

                void move_past_single_edge_node();

                void assign_end_position();

                void choose_allele();

            private:
                covG_ptr cur_Node;
                std::size_t bases_remaining;
                VariantSitePath traversed_loci;
                uint32_t traversed_index;
                bool first_node;
                node_coordinate start_pos;
                node_coordinate end_pos;
            };

            class PbCovRecorder{
            public:
                PbCovRecorder(PRG_Info& prg_info, SearchStates const& search_states,
                        std::size_t read_size);

                // Testing-related constructors
                PbCovRecorder() = default;
                PbCovRecorder(realCov_to_dummyCov existing_cov_mapping): cov_mapping(existing_cov_mapping) {}
                PbCovRecorder(PRG_Info& prg_info, std::size_t read_size) :
                        prg_info(&prg_info), read_size(read_size){}

                void process_SearchState(SearchState const& ss);
                void record_full_traversal(Traverser& t);
                void process_Node(covG_ptr cov_node, node_coordinate start_pos, node_coordinate end_pos);
                void write_coverage_from_dummy_nodes();

                realCov_to_dummyCov get_cov_mapping() const {return cov_mapping;}

            private:
                realCov_to_dummyCov cov_mapping;
                PRG_Info* prg_info; // non-const: we will add coverage in nodes of its coverage_Graph.
                std::size_t read_size;
            };
        }
    }
}
#endif //GRAMTOOLS_ALLELE_BASE_HPP
