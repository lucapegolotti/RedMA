add_executable(BypassExample BypassExample.cpp)
target_link_libraries(BypassExample ${LINKED_LIBRARIES})
target_include_directories(BypassExample PRIVATE ${HEADERS_DIRS})

add_executable(Bypass_snapshots Bypass_snapshots.cpp)
target_link_libraries(Bypass_snapshots ${LINKED_LIBRARIES})
target_include_directories(Bypass_snapshots PRIVATE ${HEADERS_DIRS})

add_executable(Bypass_generatematrices Bypass_generatematrices.cpp)
target_link_libraries(Bypass_generatematrices ${LINKED_LIBRARIES})
target_include_directories(Bypass_generatematrices PRIVATE ${HEADERS_DIRS})

add_custom_command(
    TARGET BypassExample PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/datafiles
    ${CMAKE_CURRENT_BINARY_DIR}/datafiles
)

add_custom_command(
    TARGET BypassExample PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/geometries
    ${CMAKE_CURRENT_BINARY_DIR}/geometries
)
