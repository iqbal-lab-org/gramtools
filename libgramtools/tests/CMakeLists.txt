set(INCLUDE
        ../include
        ../submods
        test_resources
        )

find_package(GTest REQUIRED)

file(GLOB_RECURSE SOURCES *.cpp)

add_executable(test_main
        main.cpp
        ${SOURCES}
        ${PROJECT_SOURCE_DIR}/libgramtools/submods/submod_resources.cpp
        ${COMMON_SOURCES} )

target_link_libraries(test_main
        gramtools
        gtest::gtest
        -lpthread
        -lm)
target_include_directories(test_main PUBLIC
        ${INCLUDE}
        )
set_target_properties(test_main
        PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON)
add_custom_command(TARGET test_main POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/test_main
        ${PROJECT_SOURCE_DIR}/libgramtools/tests/test_main.bin)


# Run test suite by issuing `make test`, or `ctest -VV`
add_test(NAME run_test_main COMMAND test_main)
