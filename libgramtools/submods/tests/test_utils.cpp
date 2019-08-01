#include "gtest/gtest.h"
#include "common/utils.hpp"


TEST(ReverseComplementRead, GivenRead_ReverseComplementReadReturned) {
    gram::Pattern read = {1, 2, 1, 3, 4};
    auto result = gram::reverse_complement_read(read);
    gram::Pattern expected = {1, 2, 4, 3, 4};
    EXPECT_EQ(result, expected);
}
