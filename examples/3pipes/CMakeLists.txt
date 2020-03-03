add_executable(3pipes_snapshots 3pipes_snapshots.cpp)
target_link_libraries(3pipes_snapshots ${LINKED_LIBRARIES})
target_include_directories(3pipes_snapshots PRIVATE ${HEADERS_DIRS})

add_executable(3pipes_generatebasis 3pipes_generatebasis.cpp)
target_link_libraries(3pipes_generatebasis ${LINKED_LIBRARIES})
target_include_directories(3pipes_generatebasis PRIVATE ${HEADERS_DIRS})

add_executable(3pipes_mdeim 3pipes_mdeim.cpp)
target_link_libraries(3pipes_mdeim ${LINKED_LIBRARIES})
target_include_directories(3pipes_mdeim PRIVATE ${HEADERS_DIRS})

add_executable(3pipes_checkmdeim 3pipes_checkmdeim.cpp)
target_link_libraries(3pipes_checkmdeim ${LINKED_LIBRARIES})
target_include_directories(3pipes_checkmdeim PRIVATE ${HEADERS_DIRS})

add_executable(3pipes_online 3pipes_online.cpp)
target_link_libraries(3pipes_online ${LINKED_LIBRARIES})
target_include_directories(3pipes_online PRIVATE ${HEADERS_DIRS})

add_custom_command(
    TARGET 3pipes_snapshots PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/datafiles
    ${CMAKE_CURRENT_BINARY_DIR}/datafiles
)

add_custom_command(
    TARGET 3pipes_snapshots PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/geometries
    ${CMAKE_CURRENT_BINARY_DIR}/geometries
)

add_custom_command(
    TARGET 3pipes_snapshots PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/scripts
    ${CMAKE_CURRENT_BINARY_DIR}/scripts
)
