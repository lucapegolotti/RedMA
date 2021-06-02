add_executable(NavierStokesNicholas NavierStokesNicholas.cpp)
target_link_libraries(NavierStokesNicholas ${LINKED_LIBRARIES})
target_include_directories(NavierStokesNicholas PRIVATE ${HEADERS_DIRS})

add_executable(NavierStokesNicholas_gen NavierStokesNicholas_gen.cpp)
target_link_libraries(NavierStokesNicholas_gen ${LINKED_LIBRARIES})
target_include_directories(NavierStokesNicholas_gen PRIVATE ${HEADERS_DIRS})

add_executable(NavierStokesNicholasImport NavierStokesNicholasImport.cpp)
target_link_libraries(NavierStokesNicholasImport ${LINKED_LIBRARIES})
target_include_directories(NavierStokesNicholasImport PRIVATE ${HEADERS_DIRS})

add_custom_command(
    TARGET NavierStokesNicholas PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/datafiles
    ${CMAKE_CURRENT_BINARY_DIR}/datafiles
)

add_custom_command(
    TARGET NavierStokesNicholas PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/geometries
    ${CMAKE_CURRENT_BINARY_DIR}/geometries
)

add_custom_command(
    TARGET NavierStokesNicholas PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/scripts
    ${CMAKE_CURRENT_BINARY_DIR}/scripts
)
