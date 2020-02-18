/** @file
 * Major role is to expose functions to select mapping instances of a read according to their equivalence classes,
 * and call coverage recording functions on selection.
 * Also has functions to set up the coverage recording structures other than coverage_Graph
 */

#ifndef GRAMTOOLS_COVERAGE_COMMON_HPP
#define GRAMTOOLS_COVERAGE_COMMON_HPP

#include "genotype/quasimap/search/types.hpp"
#include "genotype/parameters.hpp"
#include "genotype/quasimap/coverage/types.hpp"

namespace gram {

/**
 * Each type of coverage operation (record, generate, dump) operates on each level of coverage information.
 */
namespace coverage::record {
    /**
     * Selects read mappings and records all coverage information.
     * @see selection()
     */
    void search_states(Coverage &coverage,
                       const SearchStates &search_states,
                       const uint64_t &read_length,
                       const PRG_Info &prg_info,
                       const uint32_t &random_seed = 0);
}

namespace coverage::generate {
    /**
     * Calls the routines for building empty structures to record different types of coverage information.
     */
    Coverage empty_structure(const PRG_Info &prg_info);
}

namespace coverage::dump {
    /**
     * Write coverage information to disk.
     */
    void all(const Coverage &coverage,
             const GenotypeParams &parameters);
}


using SitePath = std::set<Marker>;
class RandomGenerator;
using rand_ptr = RandomGenerator *const;

/**
 * A set of site marker IDs signalling non-nested bubbles. One set defines an equivalence class.
 */
using level0_Sites = std::set<Marker>;
using uniqueLoci = std::set<VariantLocus>;

using info_ptr = PRG_Info const* const;

/**
 * Finds the set of (nested) Loci supported by a `SearchState`.
 *
 * Key attributes are:
 *  - base_sites: a `level0_Sites`. Each distinct `level0_Sites` defines an equivalence class.
 *  - unique_loci: this is a set of `VariantLocus` that the processed `SearchState` is compatible with
 *  (data struct: `uniqueLoci`).
 */
class LocusFinder{
public:

    LocusFinder() : search_state(), prg_info(nullptr) {};

    LocusFinder(SearchState const search_state, info_ptr prg_info);

    /** Sanity check: are all variant site markers in the `SearchState` different? */
    void check_site_uniqueness(SearchState const& search_state);
    void check_site_uniqueness(){check_site_uniqueness(this->search_state);}

    /**
     * Takes a `VariantLocus` and registers it as well as all sites it is nested
     * within, up to a level 0 site.
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

    level0_Sites base_sites; /**< Form the basis for `SearchState` selection */
    SitePath used_sites; /**< For remembering which sites have already been processed */
    uniqueLoci unique_loci; /**< For grouped allele counts coverage recording */
private:
    SearchState const search_state;
    info_ptr prg_info;
};

/**
 * Models an equivalence class: a list of `SearchState`s that are all compatible with the same level 0 sites.
 * The second member, `uniqueLoci`, is the set of all `VariantLocus` that the `SearchStates` are compatible with.
 */
using traversal_info = std::pair<SearchStates, uniqueLoci>;
/**
 * Models a set of equivalence classes: each `level0_Sites` is a set of site markers at level 0, ie non-nested bubbles.
 * This data structure is the basis for:
 *  -Dispatching `SearchState`s into their equivalence class
 *  -Random selection of one equivalence class.
 */
using uniqueSitePaths = std::map<level0_Sites, traversal_info>;

struct SelectedMapping{
    SearchStates navigational_search_states; /**< Use: recording per base coverage*/
    uniqueLoci equivalence_class_loci; /**< Use: recording grouped allele count and allele sum coverage*/
};

/**
 * Takes a set of `SearchState`s, dispatches them into equivalence classes, and randomly selects
 * equivalent mapping instances of the read.
 *
 * The basis for selection is the set of `level0_Sites` in `usps`.
 */
class MappingInstanceSelector{
public:
    uniqueSitePaths usps; /**< Key dispatching and selection object.*/

    // Constructor
    MappingInstanceSelector(SearchStates const search_states, info_ptr prg_info, rand_ptr rand_generator);

    // Constructors for testing
    MappingInstanceSelector() : prg_info(nullptr), rand_generator(nullptr){}
    MappingInstanceSelector(info_ptr prg_info) : prg_info(prg_info), rand_generator(nullptr){}
    MappingInstanceSelector(info_ptr prg_info, rand_ptr rand_g) : prg_info(prg_info), rand_generator(rand_g){}

    void process_searchstates(SearchStates const& all_ss);
    void set_searchstates(SearchStates const& ss) {input_search_states = ss;}

    /**
     * Dispatches a `SearchState` into `usps` using `LocusFinder`.
     */
    void add_searchstate(SearchState const& ss);

    uint32_t count_nonvar_search_states(SearchStates const& search_states);

    /**
     * Selects from the set of mapping instances of a read in the PRG.
     * Can select either:
     *  - A non-variant mapping instance: no coverage is recorded
     *  - An equivalence class of `SearchState`s: coverage is recorded for those `SearchState`s.
     */
    int32_t random_select_entry();
    void apply_selection(int32_t selected_index);
    SelectedMapping get_selection(){return selected;}
private:
    SearchStates input_search_states;
    SelectedMapping selected; /**< stores the choice made*/
    info_ptr prg_info;
    rand_ptr rand_generator;
};
}

#endif //GRAMTOOLS_COVERAGE_COMMON_HPP
