set(HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/FunctionFunctor.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SystemSolver.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SystemSolver_imp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LinearSystemSolver.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LinearSystemSolver_imp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LinearOperatorEp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/InverseOperatorEp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PreconditionerOperatorEp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SaddlePointPreconditionerEp.hpp
    CACHE INTERNAL ""
)

set(SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/LinearSystemSolver_dense.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LinearSystemSolver_sparse.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LinearOperatorEp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/InverseOperatorEp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PreconditionerOperatorEp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SaddlePointPreconditionerEp.cpp
    CACHE INTERNAL ""
)
