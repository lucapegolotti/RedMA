add_executable(TestTube TestTube.cpp ${SOURCES})
target_include_directories(TestTube PRIVATE ${HEADERS_DIRS})

configure_file(test_all.sh ${CMAKE_CURRENT_BINARY_DIR}/test_all.sh COPYONLY)
