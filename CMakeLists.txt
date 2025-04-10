cmake_minimum_required(VERSION 3.15)
project(gramtools CXX)

# Configure CCache if available - speeds up compilation of same objects with same parameters
find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

# Allows clangd code analysis in vim
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) 

# Allows `include()` of shipped modules
list(APPEND CMAKE_MODULE_PATH
     ${PROJECT_SOURCE_DIR}/libgramtools/cmake_modules)

enable_testing()
add_subdirectory(libgramtools)
