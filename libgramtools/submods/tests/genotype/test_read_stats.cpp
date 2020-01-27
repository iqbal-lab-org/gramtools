#include "gtest/gtest.h"
#include "genotype/read_stats.hpp"

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

