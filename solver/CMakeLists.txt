add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/array)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/assemblers)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/basis_functions)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/coupling)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/linear_solve)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/time_marching_algorithms)

set(HEADERS_DIRS
    ${HEADERS_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/assemblers
    ${CMAKE_CURRENT_SOURCE_DIR}/basis_functions
    ${CMAKE_CURRENT_SOURCE_DIR}/coupling
    ${CMAKE_CURRENT_SOURCE_DIR}/linear_solve
    ${CMAKE_CURRENT_SOURCE_DIR}/time_marching_algorithms
    PARENT_SCOPE
)
