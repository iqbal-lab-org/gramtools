#include "gtest/gtest.h"
#include "../../include/hello.h"
#include <sstream>

TEST(TestHello, hello) 
{
  std::ostringstream output;
  helloworld( output );
  EXPECT_STREQ("Hello World",output.str().c_str());
    }

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
