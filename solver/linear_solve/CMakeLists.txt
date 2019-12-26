set(HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/GlobalSolverOperator.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GlobalIdentityOperator.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GlobalSIMPLEOperator.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GlobalSIMPLEOperatorPseudoFSI.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SteadySolverOperator.hpp
    CACHE INTERNAL ""
)

set(SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/GlobalSolverOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GlobalIdentityOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GlobalSIMPLEOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GlobalSIMPLEOperatorPseudoFSI.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SteadySolverOperator.cpp
    CACHE INTERNAL ""
)

add_library(linear_solve SHARED
            ${CMAKE_CURRENT_SOURCE_DIR}/GlobalSolverOperator.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/GlobalSolverOperator.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/GlobalIdentityOperator.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/GlobalIdentityOperator.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/GlobalSIMPLEOperator.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/GlobalSIMPLEOperator.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/GlobalSIMPLEOperatorPseudoFSI.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/GlobalSIMPLEOperatorPseudoFSI.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/SteadySolverOperator.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/SteadySolverOperator.cpp)
