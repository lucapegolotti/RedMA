add_executable(AortaMeshing AortaMeshing.cpp)
target_link_libraries(AortaMeshing ${LINKED_LIBRARIES})
target_include_directories(AortaMeshing PRIVATE ${HEADERS_DIRS})

add_custom_command(
        TARGET AortaMeshing POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/SolverParamList.xml
                ${CMAKE_CURRENT_BINARY_DIR}/SolverParamList.xml)

add_custom_command(
        TARGET AortaMeshing PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${CMAKE_CURRENT_SOURCE_DIR}/datafiles
                ${CMAKE_CURRENT_BINARY_DIR}/datafiles)
