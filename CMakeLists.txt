cmake_minimum_required(VERSION 2.8)
project(gramtools CXX)

enable_testing()

set(SOURCE ${PROJECT_SOURCE_DIR}/src)
set(INCLUDE ${PROJECT_SOURCE_DIR}/include $ENV{HOME}/include /apps/well/htslib/1.4.1/include)
set(LIB_DIRS
        ${PROJECT_SOURCE_DIR}/lib
        $ENV{LD_LIBRARY_PATH}
        /usr/local/lib
	/apps/well/htslib/1.4.1/lib/
        $ENV{HOME}/lib)

set(CMAKE_CXX_FLAGS "-O3 -std=c++11 -Werror -Wfatal-errors -march=native -msse4.2")

SET(BOOST_ROOT /apps/well/boost/1.54.0/)
set(LIBS hts sdsl divsufsort divsufsort64 pthread z)
set(LIBS_PATHS "")
foreach(LIB ${LIBS})
    find_library(PATH_${LIB} ${LIB} HINTS ${LIB_DIRS})
    list(APPEND LIBS_PATHS ${PATH_${LIB}})
endforeach(LIB)

find_package(Boost COMPONENTS program_options timer system REQUIRED)

add_library(gram STATIC
    ${SOURCE}/map.cpp
    ${SOURCE}/kmers.cpp
    ${SOURCE}/bidir_search_bwd.cpp
    ${SOURCE}/bidir_search_fwd.cpp
    ${SOURCE}/get_location.cpp
    ${SOURCE}/skip.cpp
    ${SOURCE}/bidir_search.cpp
    ${SOURCE}/csa_construction.cpp
    ${SOURCE}/masks.cpp
    ${SOURCE}/precalc_kmer_matches.cpp
    ${SOURCE}/process_prg.cpp)
target_include_directories(gram PUBLIC
    ${SOURCE}
    ${INCLUDE}
    ${Boost_INCLUDE_DIRS})
target_link_libraries(gram LINK_PUBLIC ${LIBS_PATHS} ${Boost_LIBRARIES})

add_executable(gramtools ${SOURCE}/main.cpp)
set_target_properties(gramtools
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
target_include_directories(gramtools PUBLIC
    ${INCLUDE}
    ${Boost_INCLUDE_DIRS})
target_link_libraries(gramtools LINK_PUBLIC
    gram ${LIBS_PATHS} ${Boost_LIBRARIES})

add_subdirectory(test)
add_test(test_main test_main)
