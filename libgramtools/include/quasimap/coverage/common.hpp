#include "quasimap/search_types.hpp"
#include "quasimap/coverage/types.hpp"


#ifndef GRAMTOOLS_COVERAGE_COMMON_HPP
#define GRAMTOOLS_COVERAGE_COMMON_HPP

namespace gram {

    /**
     * Each type of coverage operation (record, generate, dump) operates on each level of coverage information.
     */
    namespace coverage {
        namespace record {
            /**
             * Selects read mappings and records coverage information.
             * @see selection()
             */
            void search_states(Coverage &coverage,
                               const SearchStates &search_states,
                               const uint64_t &read_length,
                               const PRG_Info &prg_info,
                               const uint32_t &random_seed = 0);
        }

        namespace generate {
            /**
             * Calls the routines for building empty structures to record different types of coverage information.
             */
            Coverage empty_structure(const PRG_Info &prg_info);
        }

        namespace dump {
            /**
             * Write coverage information to disk.
             */
            void all(const Coverage &coverage,
                     const Parameters &parameters);
        }
    }


    using SitePath = std::set<Marker>;
    using uniqueLoci = std::set<VariantLocus>;

    struct coverage_struct{
        SearchState search_state;
        uniqueLoci all_unique_loci;
    };
    using new_uniqueSitePaths = std::map<SitePath, std::vector<coverage_struct>>;
    using uniqueSitePaths = std::map<SitePath, SearchStates>;

    using info_ptr = PRG_Info const* const;

    /**
     * Class whose purpose it is to find the set of (nested) Loci supported by a `SearchState`
     */
    class LocusFinder{
    public:

        LocusFinder() : search_state(), prg_info(nullptr) {};

        LocusFinder(SearchState const search_state, info_ptr prg_info);

        /**
         * Takes a (potentially nested) `VariantLocus` and registers it as well as all sites it is nested
         * within, up to a base site.
         */
        void assign_nested_locus(VariantLocus const& var_loc, info_ptr info_ptr);

        /**
         * This function works on the premise that all `VariantLocus` in the `traversing_path`
         * are in the same nested bubble.
         */
        void assign_traversing_loci(SearchState const& search_state, info_ptr prg_info);
        void assign_traversing_loci(){ assign_traversing_loci(this->search_state, this->prg_info);}

        void assign_traversed_loci(SearchState const& search_state, info_ptr prg_info);
        void assign_traversed_loci(){assign_traversed_loci(this->search_state, this->prg_info);}

        SitePath base_sites; /**< 'Level 0' nesting sites; they form the basis for `SearchState` selection */
        SitePath used_sites; /**< For remembering which sites have already been processed */
        uniqueLoci unique_loci; /**< For grouped allele counts coverage recording */
    private:
        SearchState const search_state;
        info_ptr prg_info;
    };

    class MappingInstanceSelector{
    public:
        SearchStates navigational_search_states; /**< for recording per base coverage*/
        uniqueLoci equivalence_class_loci; /**< for recording grouped allele count coverage*/

        MappingInstanceSelector() : input_search_states(), usps(), prg_info(nullptr){};

        MappingInstanceSelector(SearchStates const search_states, info_ptr prg_info)
            : input_search_states(search_states), usps(), prg_info(prg_info)
            {};

        void add_searchstate(SearchState const& ss);
        // uniqueSitePaths get_unique_site_paths(const SearchStates &search_states){};
    private:
        SearchStates const input_search_states;
        info_ptr prg_info;
        new_uniqueSitePaths usps;
    };

    SitePath get_path_sites(const SearchState &search_state);
    uniqueSitePaths get_unique_site_paths(const SearchStates &search_states);

    bool check_allele_encapsulated(const SearchState &search_state,
                                   const uint64_t &read_length,
                                   const PRG_Info &prg_info);

    bool multiple_allele_encapsulated(const SearchState &search_state,
                                      const uint64_t &read_length,
                                      const PRG_Info &prg_info);

    uint64_t random_int_inclusive(const uint64_t &min,
                                  const uint64_t &max,
                                  const uint64_t &random_seed);

    uint32_t count_nonvariant_search_states(const SearchStates &search_states);



}

#endif //GRAMTOOLS_COVERAGE_COMMON_HPP
