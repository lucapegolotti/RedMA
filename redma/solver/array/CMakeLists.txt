set(HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/aMatrix.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockMatrix.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockMatrix_imp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockVector.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockVector_imp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockMaps_imp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockMaps.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/MatrixEp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/VectorEp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DenseVector.hpp
    CACHE INTERNAL ""
)

set(SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/MatrixEp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/VectorEp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockMatrix_dense.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockMatrix_sparse.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Double.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockMaps.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DenseVector.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DenseMatrix.cpp
    CACHE INTERNAL ""
)
