set(HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/BasisFunctionFunctor.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FourierBasisFunction.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ZernikeBasisFunction.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ChebyshevBasisFunction.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DummyBasisFunction.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BasisFunctionFactory.hpp
    CACHE INTERNAL ""
)

set(SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/BasisFunctionFunctor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FourierBasisFunction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ZernikeBasisFunction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ChebyshevBasisFunction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DummyBasisFunction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BasisFunctionFactory.cpp
    CACHE INTERNAL ""
)

add_library(basis_functions SHARED
            ${CMAKE_CURRENT_SOURCE_DIR}/BasisFunctionFunctor.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/BasisFunctionFunctor.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/FourierBasisFunction.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/FourierBasisFunction.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/ZernikeBasisFunction.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/ZernikeBasisFunction.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/ChebyshevBasisFunction.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/ChebyshevBasisFunction.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/DummyBasisFunction.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/DummyBasisFunction.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/BasisFunctionFactory.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/BasisFunctionFactory.cpp
)
