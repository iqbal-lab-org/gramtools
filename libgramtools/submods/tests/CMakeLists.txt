include(ExternalProject)
ExternalProject_Add(gtest
        URL https://github.com/google/googletest/archive/release-1.8.0.zip
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gtest
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

set(INCLUDE
        ${PROJECT_SOURCE_DIR}/libgramtools/include
        ${PROJECT_SOURCE_DIR}/libgramtools/submods)

file(GLOB_RECURSE SOURCES ${PROJECT_SOURCE_DIR}/libgramtools/submods/tests/*.cpp)
file(GLOB COMMON_SOURCES ${PROJECT_SOURCE_DIR}/libgramtools/submods/src_common/*.cpp)

add_executable(test_main
        main.cpp

        ${SOURCES}
        ${COMMON_SOURCES} )

target_link_libraries(test_main
        gramtools
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
        ${PROJECT_SOURCE_DIR}/cmake-build-debug/libgramtools/submods/tests/test_main
        ${PROJECT_SOURCE_DIR}/libgramtools/submods/tests)

add_test(NAME test_main COMMAND test_main)
