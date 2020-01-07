add_executable(ShowcaseDeformation ShowcaseDeformation.cpp)
target_link_libraries(ShowcaseDeformation ${LINKED_LIBRARIES})
target_include_directories(ShowcaseDeformation PRIVATE ${HEADERS_DIRS})

add_custom_command(
        TARGET ShowcaseDeformation POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/SolverParamList.xml
                ${CMAKE_CURRENT_BINARY_DIR}/SolverParamList.xml)
