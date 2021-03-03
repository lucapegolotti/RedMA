add_executable(RBPipelineSingleblock_snapshots RBPipelineSingleblock_snapshots.cpp)
target_link_libraries(RBPipelineSingleblock_snapshots ${LINKED_LIBRARIES})
target_include_directories(RBPipelineSingleblock_snapshots PRIVATE ${HEADERS_DIRS})

add_executable(RBPipelineSingleblock_generatematrices RBPipelineSingleblock_generatematrices.cpp)
target_link_libraries(RBPipelineSingleblock_generatematrices ${LINKED_LIBRARIES})
target_include_directories(RBPipelineSingleblock_generatematrices PRIVATE ${HEADERS_DIRS})

add_custom_command(
    TARGET RBPipelineSingleblock_snapshots PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/datafiles
    ${CMAKE_CURRENT_BINARY_DIR}/datafiles
)

add_custom_command(
    TARGET RBPipelineSingleblock_snapshots PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/geometries
    ${CMAKE_CURRENT_BINARY_DIR}/geometries
)

add_custom_command(
    TARGET RBPipelineSingleblock_snapshots PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/scripts
    ${CMAKE_CURRENT_BINARY_DIR}/scripts
)
