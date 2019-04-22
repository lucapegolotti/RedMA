add_executable(ExportImportXML ExportImportXML.cpp ${SOURCES})
target_link_libraries(ExportImportXML ${LINKED_LIBRARIES})
target_include_directories(ExportImportXML PRIVATE ${HEADERS_DIRS})
