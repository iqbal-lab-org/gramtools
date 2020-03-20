#include <sstream>
#include "gtest/gtest.h"

#include "genotype/infer/output_specs/segment_tracker.hpp"

using namespace gram::genotype;

class SegmentTrackerTest: public ::testing::Test{
protected:
    void SetUp(){
        std::stringstream ss1{
            "chr1\t2200\n"
            "chr2\t400\n"
        };
        std::stringstream ss2{""};
        tracker_withCoords = segmentTracker(ss1);
        tracker_noCoords = segmentTracker(ss2);
    }
    segmentTracker tracker_withCoords, tracker_noCoords;
};

TEST_F(SegmentTrackerTest, GivenNoCoordsTracker_ReturnDefaultID){
    auto result = tracker_noCoords.get_ID(1000);
    EXPECT_EQ(result, "gramtools_prg");
    result = tracker_noCoords.get_ID(40000);
    EXPECT_EQ(result, "gramtools_prg");
}

TEST_F(SegmentTrackerTest, GivenCoordsTrackerExceedBoundary_Fails){
   EXPECT_DEATH(tracker_withCoords.get_ID(40000), "");
}

TEST_F(SegmentTrackerTest, GivenCoordsTrackerPreviousSegment_Fails) {
   // Query is processed 0-based, so asking for segment boundary places us in next segment
   auto result = tracker_withCoords.get_ID(2200);
   EXPECT_EQ("chr2", result);
   // Cannot query 'backwards'
   EXPECT_DEATH(tracker_withCoords.get_ID(200), "");
}

TEST_F(SegmentTrackerTest, GivenCoordsTrackerValidQueries_GivesValidResults) {
   EXPECT_EQ(2599, tracker_withCoords.global_edge());
    EXPECT_EQ(2199, tracker_withCoords.edge());

    auto ID = tracker_withCoords.get_ID(400);
    EXPECT_EQ("chr1", ID);
    ID = tracker_withCoords.get_ID(2500);
    EXPECT_EQ("chr2", ID);
    EXPECT_EQ(2599, tracker_withCoords.edge());
}

TEST_F(SegmentTrackerTest, GivenCoordsTrackerReset_CanQueryAgain) {
    auto ID = tracker_withCoords.get_ID(2500);
    tracker_withCoords.reset();
    ID = tracker_withCoords.get_ID(100);
    EXPECT_EQ("chr1", ID);
}
