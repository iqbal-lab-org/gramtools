include(ExternalProject)
ExternalProject_Add(
        gtest
        GIT_REPOSITORY https://github.com/google/googletest/
        GIT_TAG master
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gtest
        # Disable install step
        INSTALL_COMMAND "")

# Get GTest source and binary directories from CMake project
ExternalProject_Get_Property(gtest source_dir binary_dir)

# Create a libgtest target to be used as a dependency by test programs
add_library(libgtest IMPORTED STATIC GLOBAL)
add_dependencies(libgtest gtest)

# Set libgtest properties
set_target_properties(libgtest PROPERTIES
        "IMPORTED_LOCATION" "${binary_dir}/googlemock/gtest/libgtest.a"
        "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}")

# Create a libgmock target to be used as a dependency by test programs
add_library(libgmock IMPORTED STATIC GLOBAL)
add_dependencies(libgmock gtest)

# Set libgmock properties
set_target_properties(libgmock PROPERTIES
        "IMPORTED_LOCATION" "${binary_dir}/googlemock/libgmock.a"
        "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}")

add_executable(test_main
        main.cpp
        test_kmers.cpp
        test_kmer_index.cpp
        test_bidir_search_bwd.cpp
        test_prg.cpp
        test_map.cpp)
target_link_libraries(test_main
        gramtools
        libgmock
        ${LIBS_PATHS}
        ${Boost_LIBRARIES})
target_include_directories(test_main PUBLIC
        ${SOURCE}
        ${INCLUDE}
        ${UTILS}
        ${source_dir}/googletest/include
        ${source_dir}/googlemock/include)
add_custom_command(TARGET test_main POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${PROJECT_SOURCE_DIR}/cmake-build-debug/libgramtools/tests/test_main
        ${PROJECT_SOURCE_DIR}/libgramtools/tests)

add_test(NAME test_main COMMAND test_main)
