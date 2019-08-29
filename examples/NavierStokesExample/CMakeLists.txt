add_executable(NavierStokesExample NavierStokesExample.cpp ${SOURCES})
target_link_libraries(NavierStokesExample ${LINKED_LIBRARIES})
target_include_directories(NavierStokesExample PRIVATE ${HEADERS_DIRS})

add_custom_command(
        TARGET NavierStokesExample POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/tree2.xml
                ${CMAKE_CURRENT_BINARY_DIR}/tree2.xml)

add_custom_command(
        TARGET NavierStokesExample POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/tree3.xml
                ${CMAKE_CURRENT_BINARY_DIR}/tree3.xml)

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

add_custom_command(
        TARGET NavierStokesExample POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/SolverParamList_rbf3d.xml
                ${CMAKE_CURRENT_BINARY_DIR}/SolverParamList_rbf3d.xml)

add_custom_command(
        TARGET NavierStokesExample POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/solversOptionsFast.xml
                ${CMAKE_CURRENT_BINARY_DIR}/solversOptionsFast.xml)

add_custom_command(
        TARGET NavierStokesExample POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/data
                ${CMAKE_CURRENT_BINARY_DIR}/data)
