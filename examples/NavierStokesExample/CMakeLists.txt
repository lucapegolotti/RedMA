add_executable(NavierStokesExample NavierStokesExample.cpp)
target_link_libraries(NavierStokesExample ${LINKED_LIBRARIES})
target_include_directories(NavierStokesExample PRIVATE ${HEADERS_DIRS})

add_custom_command(
        TARGET NavierStokesExample POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/tree.xml
                ${CMAKE_CURRENT_BINARY_DIR}/tree.xml)

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
                ${CMAKE_CURRENT_SOURCE_DIR}/tree4.xml
                ${CMAKE_CURRENT_BINARY_DIR}/tree4.xml)

add_custom_command(
        TARGET NavierStokesExample POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/tree5.xml
                ${CMAKE_CURRENT_BINARY_DIR}/tree5.xml)

add_custom_command(
        TARGET NavierStokesExample POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/tree6.xml
                ${CMAKE_CURRENT_BINARY_DIR}/tree6.xml)

add_custom_command(
        TARGET NavierStokesExample POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/tree7.xml
                ${CMAKE_CURRENT_BINARY_DIR}/tree7.xml)

add_custom_command(
        TARGET NavierStokesExample POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/tree8.xml
                ${CMAKE_CURRENT_BINARY_DIR}/tree8.xml)

add_custom_command(
        TARGET NavierStokesExample POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/tree9.xml
                ${CMAKE_CURRENT_BINARY_DIR}/tree9.xml)

add_custom_command(
        TARGET NavierStokesExample POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/tree10.xml
                ${CMAKE_CURRENT_BINARY_DIR}/tree10.xml)

add_custom_command(
        TARGET NavierStokesExample POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/tree11.xml
                ${CMAKE_CURRENT_BINARY_DIR}/tree11.xml)

add_custom_command(
        TARGET NavierStokesExample POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/tree12.xml
                ${CMAKE_CURRENT_BINARY_DIR}/tree12.xml)

add_custom_command(
        TARGET NavierStokesExample POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/tree13.xml
                ${CMAKE_CURRENT_BINARY_DIR}/tree13.xml)

add_custom_command(
        TARGET NavierStokesExample POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/tree14.xml
                ${CMAKE_CURRENT_BINARY_DIR}/tree14.xml)

add_custom_command(
        TARGET NavierStokesExample POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/tree15.xml
                ${CMAKE_CURRENT_BINARY_DIR}/tree15.xml)

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

add_custom_command(
        TARGET NavierStokesExample POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/open_all_vtks_paraview.py
                ${CMAKE_CURRENT_BINARY_DIR}/open_all_vtks_paraview.py)
