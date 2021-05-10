file(GLOB source_files "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

foreach (filename ${source_files})
    # creating test name
    string(FIND ${filename} "/" DASHPOS REVERSE)
    string(FIND ${filename} "." DOTPOS REVERSE)
    MATH(EXPR DASHPOS "${DASHPOS}+1")
    MATH(EXPR NAMELENGTH "${DOTPOS}-${DASHPOS}")
    string(SUBSTRING ${filename} ${DASHPOS} ${NAMELENGTH} TESTNAME)

    # linking and providing headers
    add_executable(Test_${TESTNAME} ${filename})
    target_link_libraries(Test_${TESTNAME} ${LINKED_LIBRARIES})
    target_include_directories(Test_${TESTNAME} PRIVATE ${HEADERS_DIRS})
    add_test(Test_${TESTNAME} Test_${TESTNAME})

    add_custom_command(
        TARGET Test_${TESTNAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/datafiles
        ${CMAKE_CURRENT_BINARY_DIR}/datafiles
    )
endforeach()
