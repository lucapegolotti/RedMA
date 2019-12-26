add_executable(ConvertSegmentation ConvertSegmentation.cpp)
target_link_libraries(ConvertSegmentation ${LINKED_LIBRARIES})
target_include_directories(ConvertSegmentation PRIVATE ${HEADERS_DIRS})

add_custom_command(
        TARGET ConvertSegmentation POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/SolverParamList.xml
                ${CMAKE_CURRENT_BINARY_DIR}/SolverParamList.xml)

add_custom_command(
        TARGET ConvertSegmentation PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${CMAKE_CURRENT_SOURCE_DIR}/datafiles
                ${CMAKE_CURRENT_BINARY_DIR}/datafiles)
