include(ExternalProject)

execute_process(COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/download)
execute_process(COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/src)

add_custom_target(py_git_version
        COMMAND make -C ${PROJECT_SOURCE_DIR}/gramtools/version)


###  Build type flags ##
# add a RELEASE_WITH_ASSERTS build type. Use via `cmake -DCMAKE_BUILD_TYPE=REL_WITH_ASSERTS`
# Note: this variable name needs to be uppercased here
set(CMAKE_CXX_FLAGS_REL_WITH_ASSERTS "${CMAKE_CXX_FLAGS_REL_WITH_ASSERTS} -O3")


######################
###  External libs ###
######################
include(${PROJECT_SOURCE_DIR}/libgramtools/lib/sdsl.cmake)
include(${PROJECT_SOURCE_DIR}/libgramtools/lib/htslib.cmake)
include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(TARGETS)

# Parallelisation
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNO_OPENMP")
else ()
    find_package(OpenMP REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

### Include headers ####
set(EXTERNAL_INCLUDE_DIR 
    ${CMAKE_CURRENT_BINARY_DIR}/include
    ) 

#### Link libraries ####
link_directories(${CMAKE_CURRENT_BINARY_DIR}/lib)

set(SDSL_LIBS
        -ldivsufsort
        -ldivsufsort64
        -lsdsl
        )

#########################
###  Gramtools sources ##
#########################
set(SOURCE
        ${PROJECT_SOURCE_DIR}/libgramtools/src)

file(GLOB_RECURSE SOURCE_FILES
        ${SOURCE}/*.cpp)

set(INCLUDE
        ${PROJECT_SOURCE_DIR}/libgramtools/include)

file(GLOB_RECURSE HEADER_FILES
        ${INCLUDE}/*.hpp)


#########################
####  libgramtools   ####
#########################
if (${CMAKE_EXE_LINKER_FLAGS} MATCHES "-static")
    set(STATIC_FLAGS "-static-libgcc -static-libstdc++ -ldl")
endif()
add_library(gramtools STATIC
        ${SOURCE_FILES}
        ${HEADER_FILES} # For IDE to recognise header files as part of a target
        )

target_include_directories(gramtools PUBLIC
        ${INCLUDE}
        ${EXTERNAL_INCLUDE_DIR}
        ${PROJECT_SOURCE_DIR}/libgramtools/lib
        )
target_link_libraries(gramtools LINK_PUBLIC
        ${SDSL_LIBS}
        ${CMAKE_CURRENT_BINARY_DIR}/lib/libhts.a
        CONAN_PKG::boost
        CONAN_PKG::nlohmann_json
        -lstdc++fs -lpthread -lrt -lm -lz
        ${STATIC_FLAGS}
		)
target_compile_features(gramtools PUBLIC cxx_std_17)
set_target_properties(gramtools
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
target_compile_options(gramtools PUBLIC -ftrapv -Wuninitialized)
add_dependencies(gramtools
        htslib
        py_git_version
        sdsl
        )
add_custom_command(TARGET gramtools POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory
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
        )
target_link_libraries(gram LINK_PUBLIC gramtools)
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

# Allows the release to not include backend tests
if (EXISTS "${PROJECT_SOURCE_DIR}/libgramtools/tests")
    # gcov/lcov test coverage. Activate with `cmake -DCODE_COVERAGE=ON` and run with `make test_coverage`
    # For explanation, see https://github.com/bilke/cmake-modules/blob/master/CodeCoverage.cmake
    option(CODE_COVERAGE "Enable coverage reporting" OFF)
    if(CODE_COVERAGE)
        set(COVERAGE_EXCLUDES
            "${CMAKE_CURRENT_BINARY_DIR}/include/*"
            "${CMAKE_CURRENT_BINARY_DIR}/tests/*"
            "${PROJECT_SOURCE_DIR}/libgramtools/tests/*"
            )
        include(CodeCoverage)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
        append_coverage_compiler_flags()
        setup_target_for_coverage_lcov(
            NAME test_coverage 
            DEPENDENCIES test_main
            EXECUTABLE test_main
            )
    endif()

    add_subdirectory(tests)
endif()

add_subdirectory(submods)
