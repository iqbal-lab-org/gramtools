cmake_minimum_required(VERSION 3.0.1)
project(gramtools CXX)

# Configure CCache if available - speeds up compilation of same objects with same parameters
find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

add_subdirectory(libgramtools)