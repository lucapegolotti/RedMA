add_executable(BendBifurcation BendBifurcation.cpp)
target_link_libraries(BendBifurcation ${LINKED_LIBRARIES})
target_include_directories(BendBifurcation PRIVATE ${HEADERS_DIRS})

add_custom_command(
        TARGET BendBifurcation POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/SolverParamList.xml
                ${CMAKE_CURRENT_BINARY_DIR}/SolverParamList.xml)
