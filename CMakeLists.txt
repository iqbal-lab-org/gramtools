cmake_minimum_required(VERSION 2.8)
project(gramtools CXX)

enable_testing()

set(SOURCE ${PROJECT_SOURCE_DIR}/src)
set(INCLUDE ${PROJECT_SOURCE_DIR}/include)
set(UTILS ${PROJECT_SOURCE_DIR}/utils)
set(LIB_DIR ${PROJECT_SOURCE_DIR}/lib)

set(CMAKE_CXX_FLAGS "-O3 -std=c++11 -Werror -Wfatal-errors")

set(LIBS hts sdsl divsufsort divsufsort64 pthread z)
set(LIBS_PATHS "")
foreach(LIB ${LIBS})
    find_library(PATH_${LIB} ${LIB} HINTS ${LIB_DIR} $ENV{LD_LIBRARY_PATH} /usr/local/lib $ENV{HOME}/lib)
    list(APPEND LIBS_PATHS ${PATH_${LIB}})
endforeach(LIB)
find_package(Boost COMPONENTS program_options REQUIRED)

message(STATUS ${LIBS_PATHS})

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
    ${SOURCE}/precalc_kmer_matches.cpp)
target_include_directories(gram PUBLIC
    ${SOURCE}
    ${INCLUDE}
    ${UTILS}
    ${Boost_INCLUDE_DIRS})
target_link_libraries(gram LINK_PUBLIC ${LIBS_PATHS} ${Boost_LIBRARIES})

add_executable(gramtools ${SOURCE}/main.cpp)
target_include_directories(gramtools PUBLIC
    ${INCLUDE}
    ${UTILS}
    ${Boost_INCLUDE_DIRS})
target_link_libraries(gramtools LINK_PUBLIC
    gram ${LIBS_PATHS} ${Boost_LIBRARIES})

add_subdirectory(test)
add_test(test_bidir_search test_bidir_search)
