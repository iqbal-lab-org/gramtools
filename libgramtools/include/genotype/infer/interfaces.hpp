/**
 * @file Common interface file for:
 * - genotyped sites
 * - genotyping models
 */
#ifndef INFER_IFC
#define INFER_IFC

#include <variant>
#include "genotype/infer/types.hpp"
#include "genotype/quasimap/coverage/types.hpp"

namespace gram::genotype::output_spec {
    struct site_entries;
    struct vcf_meta_info_line;
    using header_vec = std::vector<vcf_meta_info_line>;
}
using namespace gram::genotype::output_spec;

namespace gram::genotype::infer {

    /**
     * Used in allele extraction but also in level genotyper
     */
    template <typename T>
     std::vector<T> prepend(std::vector<T> const& original_object, T const& to_prepend){
        std::vector<T> result;
        result.reserve(original_object.size() + 1);
        result.insert(result.end(), to_prepend);
        result.insert(result.end(), original_object.begin(), original_object.end());

        return result;
     }

    struct gtype_information{
        allele_vector alleles;
        GtypedIndices genotype;
        allele_coverages allele_covs;
        std::size_t total_coverage; /**< Total coverage on this site */
        AlleleIds haplogroups;
    };

    /**
     * Genotyped site interface
     */
    class GenotypedSite {
    protected:
        gtype_information gtype_info;
        covG_ptr site_end_node;
        std::size_t num_haplogroups = 0; /**< The number of outgoing edges from the bubble start */

    public:
        GenotypedSite() = default;
        virtual ~GenotypedSite() {};

        gtype_information get_all_gtype_info() const {return this->gtype_info; }
        void populate_site(gtype_information const& gtype_info);
        GtypedIndices const get_genotype() const { return gtype_info.genotype; }
        allele_vector const get_alleles() const { return gtype_info.alleles; }
        covG_ptr const get_site_end_node() const { return site_end_node; }

        /** Whether the site is null genotyped */
        bool is_null() const {
            if (gtype_info.genotype.size() > 0 and gtype_info.genotype.at(0) == -1) return true;
            else return false;
        };

        void make_null() {
            gtype_info.genotype = GtypedIndices{-1};
            gtype_info.total_coverage = 0;
            this->null_model_specific_entries();
        }

        void set_alleles(allele_vector const& alleles){ gtype_info.alleles = alleles; };
        void set_genotype(GtypedIndices const& gtype){ gtype_info.genotype = gtype; }

        virtual site_entries get_model_specific_entries() = 0;
        virtual void null_model_specific_entries() = 0;

        std::size_t const &get_num_haplogroups() { return num_haplogroups; }
        bool const has_alleles() const { return gtype_info.alleles.size() > 0; }

        void set_site_end_node(covG_ptr const &end_node) { site_end_node = end_node; }
        void set_num_haplogroups(std::size_t const &num_haps) { num_haplogroups = num_haps; }

        /**
         * Given alleles and GT, return the alleles referred to by GT
         */
        allele_vector const get_unique_genotyped_alleles
                (allele_vector const &all_alleles, GtypedIndices const &genotype) const;
        allele_vector const get_unique_genotyped_alleles() const {
            return get_unique_genotyped_alleles(gtype_info.alleles, gtype_info.genotype);
        }

        /**
         * Produce the haplogroups that have not been genotyped, for use in nested
         * site invalidation.
         */
        AlleleIds const get_nonGenotyped_haplogroups() const;

        AlleleIds const get_all_haplogroups() const {
            assert(num_haplogroups > 0);
            AlleleIds result;
            for (std::size_t idx{0}; idx < num_haplogroups; idx++) result.push_back(idx);
            return result;
        }

        AlleleIds get_genotyped_haplogroups(allele_vector const& input_alleles, GtypedIndices const& input_gts) const;
    };


    /**
     * Genotyping model interface.
     * Each derived model implements the production of an abstract site.
     */
    class GenotypingModel {
        virtual gt_site_ptr get_site() = 0;
    };


    class Genotyper {
    protected:
        gt_sites genotyped_records;
        coverage_Graph const *cov_graph;
        SitesGroupedAlleleCounts const *gped_covs;
        child_map child_m;

        Genotyper() : cov_graph(nullptr), gped_covs(nullptr) {}
        Genotyper(gt_sites const& sites, child_map const& ch) :
                genotyped_records(sites), child_m(ch), cov_graph(nullptr), gped_covs(nullptr) {}

    public:
        virtual ~Genotyper() {};

        gt_sites const& get_genotyped_records() const {return genotyped_records;}
        auto const& get_cov_g() const {return cov_graph;}
        auto const& get_child_m() const {return child_m;}

        virtual header_vec get_model_specific_headers() = 0;
        void run_invalidation_process(gt_site_ptr const& genotyped_site, Marker const& site_ID);
        AlleleIds get_haplogroups_with_sites(Marker const& site_ID, AlleleIds candidate_haplogroups) const;
        void invalidate_if_needed(Marker const& parent_site_ID, AlleleIds haplogroups);
    };
}

#endif //INFER_IFC
