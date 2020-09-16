set(HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/FunctionFunctor.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SystemSolver.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LinearSystemSolver.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LinearOperator.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/InverseOperator.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PreconditionerOperator.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SaddlePointPreconditioner.hpp
    CACHE INTERNAL ""
)

set(SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/SystemSolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LinearOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/InverseOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LinearSystemSolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PreconditionerOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SaddlePointPreconditioner.cpp
    CACHE INTERNAL ""
)
