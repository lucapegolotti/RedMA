set(SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/AbstractAssembler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NavierStokesAssembler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BasisFunctionFunctor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FourierBasisFunction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ZernikeBasisFunction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GlobalSolverOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GlobalIdentityOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GlobalSIMPLEOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GlobalBlockMatrix.cpp
    CACHE INTERNAL ""
)
