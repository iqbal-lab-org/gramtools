cmake_minimum_required(VERSION 2.8)
project(gramtools CXX)

set(SOURCE ./src)
set(INCLUDE ./include)
set(UTILS ./utils)
set(LIB_DIR ./lib)
set(GTEST_DIR ./googletest/include)

include_directories(${SOURCE} ${INCLUDE} ${UTILS} ${GTEST_DIR})
set(CMAKE_CXX_FLAGS "-O3 -std=c++11 -Werror -Wfatal-errors")

set(LIBS sdsl hts divsufsort divsufsort64 pthread gtest z)
set(LIBS_PATHS "")
foreach(LIB ${LIBS})
    find_library(PATH_${LIB} ${LIB} HINTS ${LIB_DIR})
    list(APPEND LIBS_PATHS ${PATH_${LIB}})
endforeach(LIB)

find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_subdirectory(test)

add_executable(
        gramtools
        ${SOURCE}/main.cpp
        ${SOURCE}/map.cpp
        ${SOURCE}/kmers.cpp
        ${SOURCE}/bidir_search_bwd.cpp
        ${SOURCE}/bidir_search_fwd.cpp
        ${SOURCE}/get_location.cpp
        ${SOURCE}/skip.cpp
        ${SOURCE}/bidir_search.cpp
        ${SOURCE}/csa_construction.cpp
        ${SOURCE}/masks.cpp
        ${SOURCE}/precalc_kmer_matches.cpp)

target_link_libraries(gramtools ${LIBS_PATHS} ${Boost_LIBRARIES})
