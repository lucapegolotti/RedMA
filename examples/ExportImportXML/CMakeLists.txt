add_executable(ExportImportXML ExportImportXML.cpp ${SOURCES})
target_link_libraries(ExportImportXML ${LINKED_LIBRARIES})
target_include_directories(ExportImportXML PRIVATE ${HEADERS_DIRS})

add_custom_command(
        TARGET ExportImportXML POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/SolverParamList.xml
                ${CMAKE_CURRENT_BINARY_DIR}/SolverParamList.xml)
