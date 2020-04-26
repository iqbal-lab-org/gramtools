include(ExternalProject)
ExternalProject_Add(gtest
        URL https://github.com/google/googletest/archive/release-1.10.0.zip
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gtest
        INSTALL_COMMAND "")

# Get GTest source and binary directories from CMake project
ExternalProject_Get_Property(gtest source_dir binary_dir)

# Create a libgtest target to be used as a dependency by test programs
add_library(libgtest IMPORTED STATIC GLOBAL)
add_dependencies(libgtest gtest)

# Set libgtest properties
set_target_properties(libgtest PROPERTIES
        "IMPORTED_LOCATION" "${binary_dir}/lib/libgtest.a"
        "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}")

# Create a libgmock target to be used as a dependency by test programs
add_library(libgmock IMPORTED STATIC GLOBAL)
add_dependencies(libgmock gtest)

# Set libgmock properties
set_target_properties(libgmock PROPERTIES
        "IMPORTED_LOCATION" "${binary_dir}/lib/libgmock.a"
        "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}")

include_directories(${CMAKE_CURRENT_BINARY_DIR}/libgramtools/include/gtest)

set(INCLUDE
        ../include
        ../submods)

file(GLOB_RECURSE SOURCES *.cpp)

add_executable(test_main
        main.cpp
        ${SOURCES}
        ${PROJECT_SOURCE_DIR}/libgramtools/submods/submod_resources.cpp
        ${COMMON_SOURCES} )

target_link_libraries(test_main
        gramtools
        libgtest
        libgmock
        -lpthread
        -lm)
target_include_directories(test_main PUBLIC
        ${INCLUDE}
        ${source_dir}/googletest/include
        ${source_dir}/googlemock/include)
set_target_properties(test_main
        PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON)
add_custom_command(TARGET test_main POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/test_main
        ${PROJECT_SOURCE_DIR}/libgramtools/tests/test_main)

add_test(NAME test_main COMMAND test_main)
