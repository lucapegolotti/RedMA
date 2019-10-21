add_executable(RossEthierSteinmanBenchmark RossEthierSteinmanBenchmark.cpp RossEthierSteinmanDec.cpp ${SOURCES})
target_link_libraries(RossEthierSteinmanBenchmark ${LINKED_LIBRARIES})
target_include_directories(RossEthierSteinmanBenchmark PRIVATE ${HEADERS_DIRS})

add_custom_command(
        TARGET RossEthierSteinmanBenchmark PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${CMAKE_CURRENT_SOURCE_DIR}/datafiles
                ${CMAKE_CURRENT_BINARY_DIR}/datafiles)

add_custom_command(
        TARGET RossEthierSteinmanBenchmark POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/SolverParamList.xml
                ${CMAKE_CURRENT_BINARY_DIR}/SolverParamList.xml)

add_custom_command(
        TARGET RossEthierSteinmanBenchmark POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/SolverParamList2.xml
                ${CMAKE_CURRENT_BINARY_DIR}/SolverParamList2.xml)

add_custom_command(
        TARGET RossEthierSteinmanBenchmark POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/solversOptionsFast.xml
                ${CMAKE_CURRENT_BINARY_DIR}/solversOptionsFast.xml)

add_custom_command(
        TARGET RossEthierSteinmanBenchmark POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/SolverParamList_rbf3d.xml
                ${CMAKE_CURRENT_BINARY_DIR}/SolverParamList_rbf3d.xml)
