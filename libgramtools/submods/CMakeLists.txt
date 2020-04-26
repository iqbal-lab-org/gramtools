set(SUBMOD_DIR
        "${PROJECT_SOURCE_DIR}/libgramtools/submods"
        )

set(INCLUDE
        "${PROJECT_SOURCE_DIR}/libgramtools/include"
        "${SUBMOD_DIR}"
        )

set(SUBMOD_RESOURCES
        "${SUBMOD_DIR}/submod_resources.cpp"
        )

# encode_prg
add_executable(encode_prg encode_prg.cpp)
target_link_libraries(encode_prg gramtools)
target_include_directories(encode_prg PUBLIC ${INCLUDE})

add_custom_command(TARGET encode_prg POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/encode_prg
        ${SUBMOD_DIR})

# print_fm_index
add_executable(print_fm_index ${SUBMOD_RESOURCES} print_fm_index.cpp)
target_link_libraries(print_fm_index gramtools)
target_include_directories(print_fm_index PUBLIC ${INCLUDE})

add_custom_command(TARGET print_fm_index POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/print_fm_index
        ${SUBMOD_DIR})
