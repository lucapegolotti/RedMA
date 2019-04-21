add_executable(ReadXML ReadXML.cpp)

add_custom_command(
        TARGET ReadXML POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/examples/data/test.xml
                ${CMAKE_CURRENT_BINARY_DIR}/data/test.xml)


add_custom_command(
        TARGET ReadXML POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/examples/data/artery1.xml
                ${CMAKE_CURRENT_BINARY_DIR}/data/artery1.xml)

add_custom_command(
        TARGET ReadXML POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/examples/data/artery2.xml
                ${CMAKE_CURRENT_BINARY_DIR}/data/artery2.xml)

add_custom_command(
        TARGET ReadXML POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/examples/data/SolverParamList.xml
                ${CMAKE_CURRENT_BINARY_DIR}/data/SolverParamList.xml)

target_link_libraries(ReadXML ${LINKED_LIBRARIES})

add_executable(ComplexTree ComplexTree.cpp ${SOURCES})
target_link_libraries(ComplexTree ${LINKED_LIBRARIES})
target_include_directories(ComplexTree PRIVATE ${HEADERS_DIRS})

add_executable(Playground Playground.cpp ${SOURCES})
target_link_libraries(Playground ${LINKED_LIBRARIES})
target_include_directories(Playground PRIVATE ${HEADERS_DIRS})
