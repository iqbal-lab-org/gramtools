include(ExternalProject)

execute_process(COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/src)

ExternalProject_Add(sdsl
        GIT_REPOSITORY https://github.com/simongog/sdsl-lite.git
        GIT_TAG master
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}
        INSTALL_COMMAND "")
add_custom_command(TARGET sdsl
        POST_BUILD
        COMMAND bash -c "${CMAKE_CURRENT_BINARY_DIR}/src/sdsl/install.sh ${CMAKE_CURRENT_BINARY_DIR}")

ExternalProject_Add(zlib
        GIT_REPOSITORY https://github.com/madler/zlib.git
        GIT_TAG master
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}
        INSTALL_COMMAND "")
add_custom_command(TARGET zlib
        POST_BUILD
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src/zlib
        COMMAND bash -c "./configure --prefix=${CMAKE_CURRENT_BINARY_DIR}"
        COMMAND bash -c "make install")

add_custom_target(htslib
        COMMAND wget https://github.com/samtools/htslib/releases/download/1.6/htslib-1.6.tar.bz2 --timestamping --directory-prefix=${CMAKE_CURRENT_BINARY_DIR}/src
        COMMAND tar -xjf ${CMAKE_CURRENT_BINARY_DIR}/src/htslib-1.6.tar.bz2 -C ${CMAKE_CURRENT_BINARY_DIR}/src)
add_custom_command(TARGET htslib
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src/htslib-1.6
        COMMAND bash -c "autoheader"
        COMMAND bash -c "autoconf"
        COMMAND bash -c "env CPPFLAGS=-I${CMAKE_CURRENT_BINARY_DIR}/include LDFLAGS=-L${CMAKE_CURRENT_BINARY_DIR}/lib ./configure --prefix=${CMAKE_CURRENT_BINARY_DIR}"
        COMMAND bash -c "make"
        COMMAND bash -c "make install")
add_dependencies(htslib zlib)

add_custom_target(boost
        COMMAND wget https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.tar.gz --timestamping --directory-prefix=${CMAKE_CURRENT_BINARY_DIR}/src
        COMMAND tar -xvzf ${CMAKE_CURRENT_BINARY_DIR}/src/boost_1_65_1.tar.gz -C ${CMAKE_CURRENT_BINARY_DIR}/src)
add_custom_command(TARGET boost
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src/boost_1_65_1
        COMMAND ./bootstrap.sh --with-libraries=program_options,timer,system,filesystem --prefix=${CMAKE_CURRENT_BINARY_DIR}
        COMMAND ./bjam install)

add_custom_target(py_git_version
        COMMAND bash -c "make -C ${PROJECT_SOURCE_DIR}/gramtools/version")

set(CMAKE_CXX_FLAGS "-std=c++14 -Werror -Wfatal-errors -march=native -lpthread -lm -llzma -msse4.2 -g -O3")

set(EXTERNAL_LIBS_DIR
        ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(EXTERNAL_INCLUDE_DIR
        ${CMAKE_CURRENT_BINARY_DIR}/include)

set(SOURCE
        ${PROJECT_SOURCE_DIR}/libgramtools/src)
set(INCLUDE
        ${PROJECT_SOURCE_DIR}/libgramtools/include)

find_package(BZip2)

include_directories()

set(EXTERN_LIBS
        boost_filesystem
        boost_program_options
        boost_timer
        boost_system
        boost_chrono
        hts
        divsufsort
        divsufsort64
        z)

link_directories(${EXTERNAL_LIBS_DIR})

# libgramtools
add_library(gramtools STATIC
        ${SOURCE}/coverage_analysis.cpp
        ${SOURCE}/utils.cpp
        ${SOURCE}/prg.cpp
        ${SOURCE}/kmers.cpp
        ${SOURCE}/kmer_index.cpp
        ${SOURCE}/search.cpp
        ${SOURCE}/masks.cpp
        ${SOURCE}/dna_ranks.cpp
        ${SOURCE}/fm_index.cpp
        ${SOURCE}/ranks.cpp

        ${INCLUDE}/coverage_analysis.hpp
        ${INCLUDE}/utils.hpp
        ${INCLUDE}/prg.hpp
        ${INCLUDE}/kmers.hpp
        ${INCLUDE}/kmer_index.hpp
        ${INCLUDE}/kmer_index_types.hpp
        ${INCLUDE}/search.hpp
        ${INCLUDE}/masks.hpp
        ${INCLUDE}/dna_ranks.hpp
        ${INCLUDE}/fm_index.hpp
        ${INCLUDE}/ranks.hpp)
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
        ${PROJECT_SOURCE_DIR}/gramtools/lib)

# gram executable
add_executable(gram
        ${SOURCE}/main.cpp
        ${SOURCE}/timer_report.cpp)
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
