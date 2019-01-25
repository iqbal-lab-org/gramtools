/**
 * @file
 * Set of routines to dump the contents of the `gram::KmerIndex` to disk.
 * The procedures dump to disk the minimal amount of information necessary in order to fully recover the `KmerIndex` for `quasimap`.
 * The procedure works as follows:
 * * Serialise all the indexed kmers in a file `kmers` as: kmer1,kmer2,kmer3...
 * * Serialise the kmer statistics in a file `kmer_stats` as: Length_SearchStates_kmer1,Length_VariantSitePath1_kmer1,Length_VariantSitePath2_kmer1,...,Length_SearchStates_kmer2,...
 * * Serialise the kmer `gram::SearchStates` in a file `sa_intervals` as: SA_left1_kmer1,SA_right1_kmer1,SA_left2_kmer1,SA_right2_kmer1,...SA_left1_kmer2,...
 * * Serialise all their `gram::variant_site_path`s in a file `paths` as: SearchState1_Marker1_kmer1,SearchState1_Allele1_kmer1,SearchState1_Marker2_kmer1,SearchState1_Allele2_kmer1,SearchState2_Marker1_Allele1_kmer1,....
 * The `kmer_stats` file holds the information required to associate each deserialised kmer in `kmers` with its `gram::SearchStates` taken out of `search_states` and populated,
 * one by one, with `gram::variant_site_path`s from `paths`.
 */
#include "kmer_index/build.hpp"


#ifndef GRAMTOOLS_KMER_INDEX_DUMP_HPP
#define GRAMTOOLS_KMER_INDEX_DUMP_HPP

namespace gram {

    /**
     * Calculates summary statistics of the indexed kmers.
     * @ret gram::KmerIndexStats
     */
    KmerIndexStats calculate_stats(const KmerIndex &kmer_index);

    /**
     * Builds a binary file of integers ranging 0-3 representing each base of each indexed kmer.
     */
    sdsl::int_vector<3> dump_kmers(const KmerIndex &kmer_index,
                                   const Parameters &parameters);

    /**
     * Builds a binary file containing statistics for each indexed kmers.
     * These are:
     * * The number of disjoint `gram::SA_Interval`s for that kmer (ie size of `gram::SearchStates`)
     * * The length of each `gram::variant_site_path` for each `gram::SearchState` of the indexed kmer.
     * @note extracting each kmer, and each kmer's `gram::SearchStates` from @param all_kmers ensures a 1-1 mapping
     * between the dumped kmers and the dumped kmer statistics.
     */
    void dump_kmers_stats(const KmerIndexStats &stats,
                          const sdsl::int_vector<3> &all_kmers,
                          const KmerIndex &kmer_index,
                          const Parameters &parameters);

    /**
     * Builds a binary file containing the start and end SA index of each `gram::SA_Interval` for all indexed kmers.
     */
    void dump_sa_intervals(const KmerIndexStats &stats,
                           const sdsl::int_vector<3> &all_kmers,
                           const KmerIndex &kmer_index,
                           const Parameters &parameters);

    /**
     * Builds a binary file containing the variant site marker(s) and variant site allele(s) traversed for each `gram::SA_interval`.
     */
    void dump_paths(const KmerIndexStats &stats,
                    const sdsl::int_vector<3> &all_kmers,
                    const KmerIndex &kmer_index,
                    const Parameters &parameters);

    namespace kmer_index {
        /**
         * Dumps to disk the indexed kmers, their `gram::SearchStates`, their `gram::VariantSitePath`s and kmer statistics.
         */
        void dump(const KmerIndex &kmer_index, const Parameters &parameters);
    }

}

#endif //GRAMTOOLS_KMER_INDEX_DUMP_HPP
