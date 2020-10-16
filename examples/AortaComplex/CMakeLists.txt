add_executable(AortaComplex AortaComplex.cpp)
target_link_libraries(AortaComplex ${LINKED_LIBRARIES})
target_include_directories(AortaComplex PRIVATE ${HEADERS_DIRS})

add_custom_command(
        TARGET AortaComplex POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/SolverParamList.xml
                ${CMAKE_CURRENT_BINARY_DIR}/SolverParamList.xml)

add_custom_command(
        TARGET AortaComplex PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${CMAKE_CURRENT_SOURCE_DIR}/datafiles
                ${CMAKE_CURRENT_BINARY_DIR}/datafiles)

add_custom_command(
    TARGET AortaComplex PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/geometries
    ${CMAKE_CURRENT_BINARY_DIR}/geometries
)
