add_executable(TestTube TestTube.cpp)
target_include_directories(TestTube PRIVATE ${HEADERS_DIRS})
target_link_libraries(TestTube ${LINKED_LIBRARIES})

add_executable(TestBifurcationSymmetric TestBifurcationSymmetric.cpp)
target_include_directories(TestBifurcationSymmetric PRIVATE ${HEADERS_DIRS})
target_link_libraries(TestBifurcationSymmetric ${LINKED_LIBRARIES})

add_executable(TestTreeStructure TestTreeStructure.cpp)
target_include_directories(TestTreeStructure PRIVATE ${HEADERS_DIRS})
target_link_libraries(TestTreeStructure ${LINKED_LIBRARIES})

add_executable(TestGeometryParser TestGeometryParser.cpp)
target_include_directories(TestGeometryParser PRIVATE ${HEADERS_DIRS})
target_link_libraries(TestGeometryParser ${LINKED_LIBRARIES})

add_custom_command(
        TARGET TestGeometryParser POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/tests/testdata/testdata1.xml
                ${CMAKE_CURRENT_BINARY_DIR}/testdata/testdata1.xml)

configure_file(test_all.sh ${CMAKE_CURRENT_BINARY_DIR}/test_all.sh COPYONLY)
