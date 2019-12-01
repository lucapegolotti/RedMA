add_executable(RBNavierStokesSingleDomain RBNavierStokesSingleDomain.cpp ${SOURCES})
target_link_libraries(RBNavierStokesSingleDomain ${LINKED_LIBRARIES})
target_include_directories(RBNavierStokesSingleDomain PRIVATE ${HEADERS_DIRS})

add_custom_command(
        TARGET RBNavierStokesSingleDomain POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/data
                ${CMAKE_CURRENT_BINARY_DIR}/data)
