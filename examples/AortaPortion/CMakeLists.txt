add_executable(AortaPortion_snapshots AortaPortion_snapshots.cpp)
target_link_libraries(AortaPortion_snapshots ${LINKED_LIBRARIES})
target_include_directories(AortaPortion_snapshots PRIVATE ${HEADERS_DIRS})

add_executable(AortaPortion_reference AortaPortion_reference.cpp)
target_link_libraries(AortaPortion_reference ${LINKED_LIBRARIES})
target_include_directories(AortaPortion_reference PRIVATE ${HEADERS_DIRS})

add_executable(AortaPortion_generatematrices AortaPortion_generatematrices.cpp)
target_link_libraries(AortaPortion_generatematrices ${LINKED_LIBRARIES})
target_include_directories(AortaPortion_generatematrices PRIVATE ${HEADERS_DIRS})

add_executable(AortaPortion_online AortaPortion_online.cpp)
target_link_libraries(AortaPortion_online ${LINKED_LIBRARIES})
target_include_directories(AortaPortion_online PRIVATE ${HEADERS_DIRS})

add_custom_command(
    TARGET AortaPortion_snapshots PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/datafiles
    ${CMAKE_CURRENT_BINARY_DIR}/datafiles
)

add_custom_command(
    TARGET AortaPortion_snapshots PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/geometries
    ${CMAKE_CURRENT_BINARY_DIR}/geometries
)

add_custom_command(
    TARGET AortaPortion_snapshots PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/scripts
    ${CMAKE_CURRENT_BINARY_DIR}/scripts
)
