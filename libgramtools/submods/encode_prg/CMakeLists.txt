set(INCLUDE
        ${PROJECT_SOURCE_DIR}/libgramtools/include
        ${PROJECT_SOURCE_DIR}/libgramtools/submods)

file(GLOB SOURCES
        "${PROJECT_SOURCE_DIR}/libgramtools/submods/src_common/*.cpp")

add_executable(encode_prg ${SOURCES} encode_prg.cpp)
target_link_libraries(encode_prg
        gramtools )
target_include_directories(encode_prg PUBLIC
        ${INCLUDE} )

add_custom_command(TARGET encode_prg POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/encode_prg
        ${PROJECT_SOURCE_DIR}/libgramtools/submods/encode_prg)
