add_subdirectory(abstract)
add_subdirectory(block)
add_subdirectory(coupling)
add_subdirectory(finite_element)
add_subdirectory(reduced_basis)
add_subdirectory(models)

set(HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/AssemblerFactory.hpp
    CACHE INTERNAL ""
)

set(SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/AssemblerFactory.cpp
    CACHE INTERNAL ""
)
