set(HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/aMatrix.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/aMatrix_imp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/aVector.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockMatrix.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockVector.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockMaps_imp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockMaps.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SparseMatrix.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DistributedVector.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DenseVector.hpp
    CACHE INTERNAL ""
)

set(SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/SparseMatrix.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DistributedVector.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockMatrix.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockVector.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Double.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockMaps.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DenseVector.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DenseMatrix.cpp
    CACHE INTERNAL ""
)
