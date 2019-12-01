add_executable(RBNavierStokesSingleDomain RBNavierStokesSingleDomain.cpp ${SOURCES})
target_link_libraries(RBNavierStokesSingleDomain ${LINKED_LIBRARIES})
target_include_directories(RBNavierStokesSingleDomain PRIVATE ${HEADERS_DIRS})

add_custom_command(
        TARGET RBNavierStokesSingleDomain PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${CMAKE_CURRENT_SOURCE_DIR}/datafiles
                ${CMAKE_CURRENT_BINARY_DIR}/datafiles)
