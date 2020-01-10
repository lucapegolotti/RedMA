set(HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockMatrix.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AbstractMatrix.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockVector.hpp
    CACHE INTERNAL ""
)

set(SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockMatrix.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockVector.cpp
    CACHE INTERNAL ""
)
