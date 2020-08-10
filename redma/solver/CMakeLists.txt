add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/system_solver)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/time_marching_algorithms)

set(HEADERS_DIRS
    ${HEADERS_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/assemblers
    ${CMAKE_CURRENT_SOURCE_DIR}/array
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary_conditions
    ${CMAKE_CURRENT_SOURCE_DIR}/coupling
    ${CMAKE_CURRENT_SOURCE_DIR}/system_solver
    ${CMAKE_CURRENT_SOURCE_DIR}/time_marching_algorithms
    ${CMAKE_CURRENT_SOURCE_DIR}/problem
    PARENT_SCOPE
)
