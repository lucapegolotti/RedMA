add_executable(SteadyNavierStokesExample SteadyNavierStokesExample.cpp)
target_link_libraries(SteadyNavierStokesExample ${LINKED_LIBRARIES})
target_include_directories(SteadyNavierStokesExample PRIVATE ${HEADERS_DIRS})

add_executable(Bypass_SteadySnapshots Bypass_SteadySnapshots.cpp)
target_link_libraries(Bypass_SteadySnapshots ${LINKED_LIBRARIES})
target_include_directories(Bypass_SteadySnapshots PRIVATE ${HEADERS_DIRS})

add_custom_command(
    TARGET SteadyNavierStokesExample PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/datafiles
    ${CMAKE_CURRENT_BINARY_DIR}/datafiles
)

add_custom_command(
    TARGET SteadyNavierStokesExample PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/geometries
    ${CMAKE_CURRENT_BINARY_DIR}/geometries
)

add_custom_command(
    TARGET SteadyNavierStokesExample PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/scripts
    ${CMAKE_CURRENT_BINARY_DIR}/scripts
)
