add_executable(NavierStokesExample NavierStokesExample.cpp)
target_link_libraries(NavierStokesExample ${LINKED_LIBRARIES})
target_include_directories(NavierStokesExample PRIVATE ${HEADERS_DIRS})

add_custom_command(
    TARGET NavierStokesExample PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/datafiles
    ${CMAKE_CURRENT_BINARY_DIR}/datafiles
)

add_custom_command(
    TARGET NavierStokesExample PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/geometries
    ${CMAKE_CURRENT_BINARY_DIR}/geometries
)

add_custom_command(
    TARGET NavierStokesExample PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/scripts
    ${CMAKE_CURRENT_BINARY_DIR}/scripts
)
