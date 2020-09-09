add_executable(AortaIliac AortaIliac.cpp)
target_link_libraries(AortaIliac ${LINKED_LIBRARIES})
target_include_directories(AortaIliac PRIVATE ${HEADERS_DIRS})

# add_executable(AortaIliac_online AortaIliac_online.cpp)
# target_link_libraries(AortaIliac_online ${LINKED_LIBRARIES})
# target_include_directories(AortaIliac_online PRIVATE ${HEADERS_DIRS})

add_custom_command(
    TARGET AortaIliac PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/datafiles
    ${CMAKE_CURRENT_BINARY_DIR}/datafiles
)

add_custom_command(
    TARGET AortaIliac PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/geometries
    ${CMAKE_CURRENT_BINARY_DIR}/geometries
)

add_custom_command(
    TARGET AortaIliac PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/scripts
    ${CMAKE_CURRENT_BINARY_DIR}/scripts
)
