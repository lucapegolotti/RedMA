set(HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/aMatrix.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockMatrix.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockMatrix_imp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockVector.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockVector_imp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/MatrixEp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/VectorEp.hpp
    CACHE INTERNAL ""
)

set(SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/MatrixEp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/VectorEp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockMatrix.cpp
    CACHE INTERNAL ""
)
