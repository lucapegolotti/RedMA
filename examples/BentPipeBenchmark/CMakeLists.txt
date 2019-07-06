add_executable(BentPipeBenchmark BentPipeBenchmark.cpp ${SOURCES})
target_link_libraries(BentPipeBenchmark ${LINKED_LIBRARIES})
target_include_directories(BentPipeBenchmark PRIVATE ${HEADERS_DIRS})

add_custom_command(
        TARGET BentPipeBenchmark POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/tree2.xml
                ${CMAKE_CURRENT_BINARY_DIR}/tree2.xml)

add_custom_command(
        TARGET BentPipeBenchmark POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/SolverParamList.xml
                ${CMAKE_CURRENT_BINARY_DIR}/SolverParamList.xml)

add_custom_command(
        TARGET BentPipeBenchmark POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/data
                ${CMAKE_CURRENT_BINARY_DIR}/data)
