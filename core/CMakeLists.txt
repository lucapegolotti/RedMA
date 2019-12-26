set(HEADERS_DIRS
    ${HEADERS_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/geometric_structure
    ${CMAKE_CURRENT_SOURCE_DIR}/solver
    ${CMAKE_CURRENT_SOURCE_DIR}/utils
    )

add_subdirectory(geometric_structure)
add_subdirectory(solver)
add_subdirectory(utils)

add_library(core SHARED
            ${HEADERS}
            ${SOURCES})

target_link_libraries(core ${LINKED_LIBRARIES})
