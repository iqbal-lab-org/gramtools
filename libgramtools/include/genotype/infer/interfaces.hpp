/**
 * @file Common interface file for:
 * - genotyped sites
 * - genotyping models
 */
#ifndef INFER_IFC
#define INFER_IFC

#include <variant>
#include "genotype/infer/types.hpp"

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

    using GtypedIndex = std::size_t; /**< The index of an allele in an allele vector */
    using GtypedIndices = std::vector<GtypedIndex>;
    using GenotypeOrNull = std::variant<GtypedIndices, bool>;
    using allele_coverages = std::vector<double>;

    class AbstractGenotypedSite;
    using gt_site = AbstractGenotypedSite;
    using gt_site_ptr = std::shared_ptr<AbstractGenotypedSite>;
    using gt_sites = std::vector<gt_site_ptr>;

    /**
     * Genotyped site interface
     */
    class AbstractGenotypedSite {
    protected:
        allele_vector alleles;
        GenotypeOrNull genotype;
        covG_ptr site_end_node;
        std::size_t num_haplogroups = 0; /**< The number of outgoing edges from the bubble start */

    public:
        virtual ~AbstractGenotypedSite() {};

        virtual GenotypeOrNull const get_genotype() const = 0;

        virtual allele_vector const get_alleles() const = 0;

        virtual covG_ptr const get_site_end_node() const = 0;

        virtual bool is_null() const = 0;

        virtual void make_null() = 0;

        std::size_t const &get_num_haplogroups() { return num_haplogroups; }
        bool const has_alleles() const { return alleles.size() > 0; }

        void set_site_end_node(covG_ptr const &end_node) { site_end_node = end_node; }
        void set_num_haplogroups(std::size_t const &num_haps) { num_haplogroups = num_haps; }


        /**
         * Given alleles and GT, return the alleles referred to by GT
         */
        allele_vector const get_unique_genotyped_alleles
                (allele_vector const &all_alleles, GenotypeOrNull const &genotype) const;

        /**
         * Return site's called alleles directly, from its own alleles
         */
        allele_vector const get_unique_genotyped_alleles() const {
            return get_unique_genotyped_alleles(alleles, genotype);
        }

        /**
         * Return site's called alleles indirectly
         * This version exists for allowing to use mocked alleles and genotypes
         */
        allele_vector const extract_unique_genotyped_alleles() const {
            auto extracted_alleles = this->get_alleles();
            auto extracted_gts = this->get_genotype();
            return get_unique_genotyped_alleles(extracted_alleles, extracted_gts);
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
    };


    /**
     * Genotyping model interface.
     * Each derived model implements the production of an abstract site.
     */
    class AbstractGenotypingModel {
        virtual gt_site_ptr get_site() = 0;
    };

}

#endif //INFER_IFC
