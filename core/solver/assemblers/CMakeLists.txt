set(HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/AbstractAssembler.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AssemblersFactory.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AssemblersFactory.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GlobalAssembler.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NavierStokesAssembler.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SUPGStabilization.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PseudoFSIAssembler.hpp
    CACHE INTERNAL ""
)

set(SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/AbstractAssembler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AssemblersFactory.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GlobalAssembler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NavierStokesAssembler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SUPGStabilization.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PseudoFSIAssembler.cpp
    CACHE INTERNAL ""
)
