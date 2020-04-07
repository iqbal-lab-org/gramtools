set(INCLUDE
        ${PROJECT_SOURCE_DIR}/libgramtools/include
        ${PROJECT_SOURCE_DIR}/libgramtools/submods)

file(GLOB SOURCES
        "${PROJECT_SOURCE_DIR}/libgramtools/submods/src_common/*.cpp")

add_executable(print_fm_index ${SOURCES} print_fm_index.cpp)
target_link_libraries(print_fm_index
        gramtools)
target_include_directories(print_fm_index PUBLIC
        ${INCLUDE} )

add_custom_command(TARGET print_fm_index POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/print_fm_index
        ${PROJECT_SOURCE_DIR}/libgramtools/submods/print_fm_index)
