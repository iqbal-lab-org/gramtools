/** @file
 *  Compute read-related statistics.
 *  For inferring a personalised reference, we need estimation of base-level error-rate, and estimation of expected coverage.
 *  This module allows for recording that, as well as some other usable metrics, such as max read length and number of
 *  sites with no coverage.
 */
#include "genotype/quasimap/coverage/types.hpp"

#ifndef GRAMTOOLS_READSTATS_HPP
#define GRAMTOOLS_READSTATS_HPP

#define NUM_READS_USED 10000

namespace gram {

    /**
     * Class for recording read-related statistics.
     */
    class ReadStats {
    public:
        // Default constructor: -1 initialisation to signal that attribute has not been computed.
        ReadStats() : mean_cov_depth(-1), no_qual_reads(-1), max_read_length(0), num_bases_processed(-1),
                      mean_pb_error(-1), variance_depth(-1), num_sites_noCov(0), num_sites_total(-1) {};

        /**
         * From a file
         */
        void compute_base_error_rate(const std::string &reads_fpath);

        /**
         * From random access memory
         */
        void compute_base_error_rate(GenomicRead_vector const& reads);

        /**
         * Compute probability of erroneous base from base Phred scores.
         */
        void process_read_perbase_error_rates(AbstractGenomicReadIterator& reads_it);


        /**
         * Compute the depth of coverage using recorded coverage of reads over variant sites after `quasimap`.
         * @param coverage gram::Coverage containing read coverage over variant sites.
         */
        void compute_coverage_depth(Coverage const& coverage, parental_map const &par_map);

        void serialise(const std::string &json_output_fpath);

        double const& get_mean_pb_error() const{ return mean_pb_error; }
        int64_t const& get_num_bases_processed() const{ return num_bases_processed; }
        std::size_t const& get_max_read_len() const{ return max_read_length; }
        int64_t const& get_num_no_qual_reads() const{ return no_qual_reads; }

        double const& get_mean_cov_depth() const{ return mean_cov_depth; }
        std::size_t const& get_num_sites_noCov() const{ return num_sites_noCov; }
        std::size_t const& get_num_sites_total() const{ return num_sites_total; }


    private:
        double mean_pb_error; // Pb sequencing error rate
        int64_t no_qual_reads;
        std::size_t max_read_length;
        int64_t num_bases_processed;

        double mean_cov_depth; // Mean of total coverages at each variant site
        double variance_depth;
        std::size_t num_sites_noCov;
        std::size_t num_sites_total;
    };

}

#endif // GRAMTOOLS_READSTATS_HPP
