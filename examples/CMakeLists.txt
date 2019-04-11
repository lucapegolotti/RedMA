add_executable(ReadXML ReadXML.cpp)

add_custom_command(
        TARGET ReadXML POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/examples/data/test.xml
                ${CMAKE_CURRENT_BINARY_DIR}/data/test.xml)

target_link_libraries(ReadXML ${LINKED_LIBRARIES})

add_executable(InstantiateTube InstantiateTube.cpp ${SOURCES})
target_include_directories(InstantiateTube PRIVATE ${HEADERS_DIRS})
target_link_libraries(InstantiateTube ${LINKED_LIBRARIES})
