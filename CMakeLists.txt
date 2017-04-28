cmake_minimum_required(VERSION 2.8)
project(gramtools CXX)

set(SOURCE ./src)
set(INCLUDE ./include)
set(UTILS ./utils)
set(LIB_DIR ./lib)
set(GTEST_DIR ./googletest/include)

include_directories(${SOURCE} ${INCLUDE} ${UTILS} ${GTEST_DIR})

set(CMAKE_CXX_FLAGS "-O3 -std=c++11 -Werror -Wfatal-errors")
add_executable(
        gramtools
        ${SOURCE}/map.cpp
        ${SOURCE}/bidir_search_bwd.cpp
        ${SOURCE}/bidir_search_fwd.cpp
        ${SOURCE}/get_location.cpp
        ${SOURCE}/precalc_gen.hpp
        ${SOURCE}/skip.cpp
        ${SOURCE}/bidir_search.cpp
        ${SOURCE}/csa_construction.cpp
        ${SOURCE}/masks.cpp
        ${SOURCE}/precalc_kmer_matches.cpp)

set(LIBS sdsl hts divsufsort divsufsort64 pthread gtest z)
SET(LIBS_PATHS "")
FOREACH(LIB ${LIBS})
    FIND_LIBRARY(PATH_${LIB} ${LIB} HINTS ${LIB_DIR})
    LIST(APPEND LIBS_PATHS ${PATH_${LIB}})
ENDFOREACH(LIB)

find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

target_link_libraries(gramtools ${LIBS_PATHS} ${Boost_LIBRARIES})
