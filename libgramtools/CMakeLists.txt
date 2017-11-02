add_custom_target(git_version
                  COMMAND bash -c
                  "make -C ${PROJECT_SOURCE_DIR}/libgramtools/include/git_version")

set(SOURCE
        ${PROJECT_SOURCE_DIR}/libgramtools/src)
set(INCLUDE
        ${PROJECT_SOURCE_DIR}/libgramtools/include
        $ENV{HOME}/include
        /apps/well/htslib/1.4.1/include)
set(LIB_DIRS
        ${PROJECT_SOURCE_DIR}/libgramtools/lib
        $ENV{LD_LIBRARY_PATH}
        /usr/local/lib
	    /apps/well/htslib/1.4.1/lib/
        $ENV{HOME}/lib)

set(CMAKE_CXX_FLAGS "-std=c++14 -Werror -Wfatal-errors -march=native -msse4.2 -g -O3")

set(BOOST_ROOT /apps/well/boost/1.54.0/)
set(LIBS
        hts
        sdsl
        divsufsort
        divsufsort64
        pthread
        z)
set(LIBS_PATHS "")
foreach(LIB ${LIBS})
    find_library(PATH_${LIB} ${LIB} HINTS ${LIB_DIRS})
    list(APPEND LIBS_PATHS ${PATH_${LIB}})
endforeach(LIB)

find_package(Boost COMPONENTS
        program_options
        timer
        system REQUIRED)

# libgramtools
add_library(gramtools STATIC
    ${SOURCE}/coverage_analysis.cpp
    ${SOURCE}/utils.cpp
    ${SOURCE}/prg.cpp
    ${SOURCE}/kmer_index.cpp
    ${SOURCE}/search.cpp
    ${SOURCE}/masks.cpp
    ${SOURCE}/fm_index.cpp
    ${SOURCE}/ranks.cpp)
target_include_directories(gramtools PUBLIC
    ${SOURCE}
    ${INCLUDE}
    ${Boost_INCLUDE_DIRS})
target_link_libraries(gramtools LINK_PUBLIC
        ${LIBS_PATHS}
        ${Boost_LIBRARIES})
set_target_properties(gramtools
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake-build-debug/lib
        LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake-build-debug/lib)
add_dependencies(gramtools git_version)
add_custom_command(TARGET gramtools POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory
            ${PROJECT_SOURCE_DIR}/gramtools/lib
        COMMAND ${CMAKE_COMMAND} -E copy
            ${PROJECT_SOURCE_DIR}/cmake-build-debug/lib/libgramtools.a
            ${PROJECT_SOURCE_DIR}/gramtools/lib)

# gram executable
add_executable(gram
        ${SOURCE}/main.cpp
        ${SOURCE}/timer_report.cpp)
target_include_directories(gram PUBLIC
    ${INCLUDE}
    ${Boost_INCLUDE_DIRS})
target_link_libraries(gram LINK_PUBLIC
        gramtools
        ${LIBS_PATHS}
        ${Boost_LIBRARIES})
set_target_properties(gram
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake-build-debug/bin)
add_custom_command(TARGET gram POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory
            ${PROJECT_SOURCE_DIR}/gramtools/bin
        COMMAND ${CMAKE_COMMAND} -E copy
            ${PROJECT_SOURCE_DIR}/cmake-build-debug/bin/gram
            ${PROJECT_SOURCE_DIR}/gramtools/bin)

# tests
enable_testing()
add_subdirectory(tests)
add_test(test_main test_main)
