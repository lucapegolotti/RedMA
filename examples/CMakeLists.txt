add_executable(ReadXML ReadXML.cpp)

add_custom_command(
        TARGET ReadXML POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/examples/data/test.xml
                ${CMAKE_CURRENT_BINARY_DIR}/data/test.xml)

target_link_libraries(ReadXML ${LINKED_LIBRARIES})
