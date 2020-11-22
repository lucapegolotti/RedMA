file(GLOB allfiles "${CMAKE_CURRENT_SOURCE_DIR}/*")

set(HEADERS_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/AtomicTest.hpp" CACHE INTERNAL "")
set(SOURCES_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/AtomicTest.cpp" CACHE INTERNAL "")

add_library(tests STATIC ${HEADERS_TESTS} ${SOURCES_TESTS})
target_link_libraries(tests ${LINKED_LIBRARIES})

set(LINKED_LIBRARIES ${LINKED_LIBRARIES} tests)

foreach(filename ${allfiles})
    if(IS_DIRECTORY ${filename})
        add_subdirectory(${filename})
    endif()
endforeach()
