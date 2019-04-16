add_executable(TestTube TestTube.cpp ${SOURCES})
target_include_directories(TestTube PRIVATE ${HEADERS_DIRS})
target_link_libraries(TestTube ${LINKED_LIBRARIES})

add_executable(TestBifurcationSymmetric TestBifurcationSymmetric.cpp ${SOURCES})
target_include_directories(TestBifurcationSymmetric PRIVATE ${HEADERS_DIRS})
target_link_libraries(TestBifurcationSymmetric ${LINKED_LIBRARIES})

add_executable(TestTreeStructure TestTreeStructure.cpp ${SOURCES})
target_include_directories(TestTreeStructure PRIVATE ${HEADERS_DIRS})
target_link_libraries(TestTreeStructure ${LINKED_LIBRARIES})

configure_file(test_all.sh ${CMAKE_CURRENT_BINARY_DIR}/test_all.sh COPYONLY)
