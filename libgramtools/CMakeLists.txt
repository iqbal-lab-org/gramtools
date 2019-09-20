include(ExternalProject)

execute_process(COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/download)
execute_process(COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/src)

include(${PROJECT_SOURCE_DIR}/libgramtools/lib/boost.cmake)
include(${PROJECT_SOURCE_DIR}/libgramtools/lib/sdsl.cmake)
include(${PROJECT_SOURCE_DIR}/libgramtools/lib/zlib.cmake)
include(${PROJECT_SOURCE_DIR}/libgramtools/lib/htslib.cmake)
add_dependencies(htslib zlib)

add_custom_target(py_git_version
        COMMAND make -C ${PROJECT_SOURCE_DIR}/gramtools/version)

set(CMAKE_CXX_FLAGS "-std=c++17 -march=native -lpthread -lrt -lm -llzma -msse4.2 -fopenmp -ftrapv -g -O3")

set(EXTERNAL_LIBS_DIR
        ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(EXTERNAL_INCLUDE_DIR
        ${CMAKE_CURRENT_BINARY_DIR}/include)

set(SOURCE
        ${PROJECT_SOURCE_DIR}/libgramtools/src)
set(INCLUDE
        ${PROJECT_SOURCE_DIR}/libgramtools/include)

find_package(BZip2)

# In development mode, makes sure can access the dependencies
include_directories(${EXTERNAL_INCLUDE_DIR})

set(EXTERN_LIBS
        boost_filesystem
        boost_program_options
        boost_timer
        boost_system
        boost_chrono
        boost_random
        boost_serialization
        hts
        divsufsort
        divsufsort64
        z)

link_directories(${EXTERNAL_LIBS_DIR})

file(GLOB_RECURSE SOURCE_FILES
        ${SOURCE}/*.cpp)

file(GLOB_RECURSE INCLUDE_FILES
        ${INCLUDE}/*.hpp)

# libgramtools
add_library(gramtools STATIC
        ${SOURCE_FILES}
        ${INCLUDE_FILES})
target_include_directories(gramtools PUBLIC
        ${INCLUDE}
        ${EXTERNAL_INCLUDE_DIR}
        ${BZIP_INCLUDE_DIRS}
        ${LZMA_INCLUDE_DIR})
target_link_libraries(gramtools LINK_PUBLIC
        ${EXTERN_LIBS}
        ${BZIP2_LIBRARIES}
        ${CMAKE_CURRENT_BINARY_DIR}/lib/libsdsl.a)
set_target_properties(gramtools
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake-build-debug/lib
        LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake-build-debug/lib)
add_dependencies(gramtools
        sdsl
        zlib
        htslib
        boost
        py_git_version)
add_custom_command(TARGET gramtools POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory
        ${PROJECT_SOURCE_DIR}/gramtools/lib
        COMMAND ${CMAKE_COMMAND} -E copy
        ${PROJECT_SOURCE_DIR}/cmake-build-debug/lib/libgramtools.a
        ${PROJECT_SOURCE_DIR}/gramtools/lib
        COMMAND ${CMAKE_COMMAND} -E copy
        ${PROJECT_SOURCE_DIR}/cmake-build-debug/libgramtools/lib/*
        ${PROJECT_SOURCE_DIR}/gramtools/lib)

# gram executable
add_executable(gram
        ${SOURCE}/main.cpp
        ${SOURCE}/common/timer_report.cpp)
add_dependencies(gram py_git_version)
target_include_directories(gram PUBLIC
        ${INCLUDE}
        ${EXTERNAL_INCLUDE_DIR}
        ${BZIP_INCLUDE_DIRS}
        ${LZMA_INCLUDE_DIR})
target_link_libraries(gram LINK_PUBLIC
        gramtools
        ${EXTERN_LIBS}
        ${BZIP2_LIBRARIES}
        ${CMAKE_CURRENT_BINARY_DIR}/lib/libsdsl.a)
set_target_properties(gram
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake-build-debug/bin
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON)
add_custom_command(TARGET gram POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory
        ${PROJECT_SOURCE_DIR}/gramtools/bin
        COMMAND ${CMAKE_COMMAND} -E copy
        ${PROJECT_SOURCE_DIR}/cmake-build-debug/bin/gram
        ${PROJECT_SOURCE_DIR}/gramtools/bin)

# tests
enable_testing()
add_subdirectory(submods/tests)
add_test(test_main test_main)

add_subdirectory(submods/print_fm_index)