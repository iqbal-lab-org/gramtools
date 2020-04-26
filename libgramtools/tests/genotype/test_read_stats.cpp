#include "gtest/gtest.h"

#include "genotype/read_stats.hpp"

#include "../test_resources/test_resources.hpp"

using namespace gram;

TEST(ReadProcessingStats, GivenTwoGenomicReads_CorrectReadProcessingStats){
    GenomicRead_vector reads{
        GenomicRead{"Read1", "AAAA", "5555"}, // 5: ASCII 53, Q-score 20 (Phred +33 scale), error prob 0.01.
        GenomicRead{"Read1", "TTTT", "5555"}
    };

    ReadStats r;
    r.compute_base_error_rate(reads);

    EXPECT_EQ(r.get_num_bases_processed(), 8);
    EXPECT_EQ(r.get_max_read_len(), 4);
    EXPECT_FLOAT_EQ(r.get_mean_pb_error(), 0.01);
}

TEST(ReadProcessingStats, GivenOneOKAndOneEmptyGenomicRead_CorrectReadProcessingStats){
    GenomicRead_vector reads{
            GenomicRead{"Read1", "AAA", "???"}, // ASCII 63, Q-score 30, error prob 0.001
            GenomicRead{"Read1", "", ""}
    };

    ReadStats r;
    r.compute_base_error_rate(reads);

    EXPECT_EQ(r.get_num_no_qual_reads(), 1);
    EXPECT_FLOAT_EQ(r.get_mean_pb_error(), 0.001);
}


/**
* NOTE: Because coverage is 'propagated upwards', we only want to record mapping statistics
* on level 1 sites (sites not nested in any others)
 */
TEST(ReadMappingStats, GivenCovStructAndParMap_CorrectMappingRelatedStats) {
    ReadStats stats;
    Coverage cov{
            {},
            SitesGroupedAlleleCounts {
                    GroupedAlleleCounts { {AlleleIds{0}, 20 }},
                    GroupedAlleleCounts { {AlleleIds{1}, 2 }},
                    GroupedAlleleCounts {}
            },
            {}
    };
    // Here we are saying that site with index 1 is nested within site with index 0
    parental_map par_map{
            {7, VariantLocus{5, 2}}
    };
    stats.compute_coverage_depth(cov, par_map);

    EXPECT_FLOAT_EQ(stats.get_mean_cov_depth(), 10); // mean of 20 and 0
    EXPECT_EQ(stats.get_num_sites_noCov(), 1); // The site with index 2 has no cov
    EXPECT_EQ(stats.get_num_sites_total(), 2); // There are two level 1 sites and one nested one.
}

/**
 * We can also test cov stats this (integration) way. It is more costly, but it relies
 * on proper mapping, coverage recording and parental map formation and thus also tests those.
 */
TEST(ReadMappingStats, GivenThreeMappedReadsNonNestedPRG_CorrectMappingRelatedStats){
    GenomicRead_vector reads{
            GenomicRead{"Read1", "AAA", "###"}, // '#' = Q-score of 2
            GenomicRead{"Read2", "AAAA", "####"},
            GenomicRead{"Read3", "GAAAA", "#####"}
    };

    prg_setup setup;
    Sequences kmers{ encode_dna_bases("AA") };
    setup.setup_numbered_prg("G5AAAA6AA6T7G8C8GGG", kmers);
    setup.quasimap_reads(reads);

    auto stats = setup.read_stats;
    // Map 3 reads to site 1, and 0 to site 2, so expect mean of (3 + 0) / 2
    EXPECT_FLOAT_EQ(stats.get_mean_cov_depth(),1.5);
    EXPECT_EQ(stats.get_num_sites_noCov(), 1);
    EXPECT_EQ(stats.get_num_sites_total(), 2);
}

TEST(ReadMappingStats, GivenTwoMappedReadsNestedPRG_CorrectMappingRelatedStats){
    GenomicRead_vector reads{
            GenomicRead{"Read1", "GGGGGCCC", "IIIIIIIII"}, // 'I' = Q-score of 40
            GenomicRead{"Read2", "GCCC", "IIII"},
    };

    prg_setup setup;
    Sequences kmers{ encode_dna_bases("CC") };
    setup.setup_bracketed_prg("G[GG[G,A]G,C]CCC", kmers);
    setup.quasimap_reads(reads);

    auto stats = setup.read_stats;
    EXPECT_FLOAT_EQ(stats.get_mean_cov_depth(), 2);
    EXPECT_EQ(stats.get_num_sites_noCov(), 0);
    EXPECT_EQ(stats.get_num_sites_total(), 1);
}
