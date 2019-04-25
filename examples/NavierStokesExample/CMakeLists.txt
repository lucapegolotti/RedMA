add_executable(NavierStokesExample NavierStokesExample.cpp ${SOURCES})
target_link_libraries(NavierStokesExample ${LINKED_LIBRARIES})
target_include_directories(NavierStokesExample PRIVATE ${HEADERS_DIRS})

add_custom_command(
        TARGET NavierStokesExample POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/tree7.xml
                ${CMAKE_CURRENT_BINARY_DIR}/tree7.xml)

add_custom_command(
        TARGET NavierStokesExample POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/SolverParamList.xml
                ${CMAKE_CURRENT_BINARY_DIR}/SolverParamList.xml)
