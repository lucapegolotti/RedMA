add_executable(ReadXML ReadXML.cpp)
target_link_libraries(ReadXML ${LINKED_LIBRARIES})

add_custom_command(
        TARGET ReadXML POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/test.xml
                ${CMAKE_CURRENT_BINARY_DIR}/test.xml)
