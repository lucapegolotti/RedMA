add_executable(AortaPortion_snapshots AortaPortion_snapshots.cpp)
target_link_libraries(AortaPortion_snapshots ${LINKED_LIBRARIES})
target_include_directories(AortaPortion_snapshots PRIVATE ${HEADERS_DIRS})

add_executable(AortaPortion_generatebasis AortaPortion_generatebasis.cpp)
target_link_libraries(AortaPortion_generatebasis ${LINKED_LIBRARIES})
target_include_directories(AortaPortion_generatebasis PRIVATE ${HEADERS_DIRS})

add_executable(AortaPortion_mdeim AortaPortion_mdeim.cpp)
target_link_libraries(AortaPortion_mdeim ${LINKED_LIBRARIES})
target_include_directories(AortaPortion_mdeim PRIVATE ${HEADERS_DIRS})

add_executable(AortaPortion_checkmdeim AortaPortion_checkmdeim.cpp)
target_link_libraries(AortaPortion_checkmdeim ${LINKED_LIBRARIES})
target_include_directories(AortaPortion_checkmdeim PRIVATE ${HEADERS_DIRS})

add_executable(AortaPortion_online AortaPortion_online.cpp)
target_link_libraries(AortaPortion_online ${LINKED_LIBRARIES})
target_include_directories(AortaPortion_online PRIVATE ${HEADERS_DIRS})

add_executable(AortaPortion_reference AortaPortion_reference.cpp)
target_link_libraries(AortaPortion_reference ${LINKED_LIBRARIES})
target_include_directories(AortaPortion_reference PRIVATE ${HEADERS_DIRS})

add_executable(AortaPortion_femvsrb AortaPortion_femvsrb.cpp)
target_link_libraries(AortaPortion_femvsrb ${LINKED_LIBRARIES})
target_include_directories(AortaPortion_femvsrb PRIVATE ${HEADERS_DIRS})

add_executable(AortaPortion_snapshotsdeim AortaPortion_snapshotsdeim.cpp)
target_link_libraries(AortaPortion_snapshotsdeim ${LINKED_LIBRARIES})
target_include_directories(AortaPortion_snapshotsdeim PRIVATE ${HEADERS_DIRS})

add_executable(AortaPortion_generatematricesforbasis AortaPortion_generatematricesforbasis.cpp)
target_link_libraries(AortaPortion_generatematricesforbasis ${LINKED_LIBRARIES})
target_include_directories(AortaPortion_generatematricesforbasis PRIVATE ${HEADERS_DIRS})

add_executable(AortaPortion_dumpsnapshots AortaPortion_dumpsnapshots.cpp)
target_link_libraries(AortaPortion_dumpsnapshots ${LINKED_LIBRARIES})
target_include_directories(AortaPortion_dumpsnapshots PRIVATE ${HEADERS_DIRS})

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
