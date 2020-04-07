include(ExternalProject)

execute_process(COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/download)
execute_process(COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/src)

add_custom_target(py_git_version
        COMMAND make -C ${PROJECT_SOURCE_DIR}/gramtools/version)


######################
###  Compile flags ##
######################
# add a RELEASE_WITH_ASSERTS build type. Use via `cmake -DCMAKE_BUILD_TYPE=REL_WITH_ASSERTS`
# Note: this variable name NEEDS TO BE UPPERCASE (but not in cmake call)
set(CMAKE_CXX_FLAGS_REL_WITH_ASSERTS "-O3")


######################
###  External libs ###
######################
include(${PROJECT_SOURCE_DIR}/libgramtools/lib/boost.cmake)
include(${PROJECT_SOURCE_DIR}/libgramtools/lib/sdsl.cmake)
include(${PROJECT_SOURCE_DIR}/libgramtools/lib/zlib.cmake)
include(${PROJECT_SOURCE_DIR}/libgramtools/lib/htslib.cmake)
include(${PROJECT_SOURCE_DIR}/libgramtools/lib/json.cmake)
add_dependencies(htslib zlib)

## Boost ##
find_package(BZip2)
find_package(OpenMP REQUIRED)

### Include headers ####
set(EXTERNAL_INCLUDE_DIR
        ${CMAKE_CURRENT_BINARY_DIR}/include)

#### Link libraries ####
link_directories(${CMAKE_CURRENT_BINARY_DIR}/lib)

set(EXTERN_LIBS
        boost_timer
        boost_chrono
        boost_program_options
        boost_serialization
        boost_system
        hts
        divsufsort
        divsufsort64
        z)

#########################
###  Gramtools sources ##
#########################
set(SOURCE
        ${PROJECT_SOURCE_DIR}/libgramtools/src)
set(INCLUDE
        ${PROJECT_SOURCE_DIR}/libgramtools/include)


file(GLOB_RECURSE SOURCE_FILES
        ${SOURCE}/*.cpp)

file(GLOB_RECURSE INCLUDE_FILES
        ${INCLUDE}/*.hpp)

#########################
####  libgramtools   ####
#########################
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
        ${CMAKE_CURRENT_BINARY_DIR}/lib/libsdsl.a
        OpenMP::OpenMP_CXX
        -lstdc++fs -lpthread -lrt -lm -llzma
		)
target_compile_features(gramtools PUBLIC cxx_std_17)
set_target_properties(gramtools
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
target_compile_options(gramtools PUBLIC -ftrapv -Wuninitialized)
add_dependencies(gramtools
        boost
        htslib
        nlohmann_json
        py_git_version
        sdsl
        zlib)
add_custom_command(TARGET gramtools POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory
        ${PROJECT_SOURCE_DIR}/gramtools/lib
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/lib/libgramtools.a
        ${PROJECT_SOURCE_DIR}/gramtools/lib
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/lib/*
        ${PROJECT_SOURCE_DIR}/gramtools/lib)

#########################
###  gram executable  ###
#########################
add_executable(gram
        ${SOURCE}/main.cpp
        ${SOURCE}/common/timer_report.cpp)
add_dependencies(gram py_git_version)
target_include_directories(gram PUBLIC
        ${INCLUDE}
        ${EXTERNAL_INCLUDE_DIR}
        ${BZIP_INCLUDE_DIRS}
        ${LZMA_INCLUDE_DIR})
target_link_libraries(gram LINK_PUBLIC gramtools )
set_target_properties(gram
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON)
add_custom_command(TARGET gram POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory
        ${PROJECT_SOURCE_DIR}/gramtools/bin
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/bin/gram
        ${PROJECT_SOURCE_DIR}/gramtools/bin)

#################
####  tests  ####
#################
enable_testing()
add_subdirectory(submods/tests)
add_test(test_main test_main)

add_subdirectory(submods/print_fm_index)
add_subdirectory(submods/encode_prg)
