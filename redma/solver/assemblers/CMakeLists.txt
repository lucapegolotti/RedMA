set(HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/AssemblerFactory.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/aAssembler.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/StokesAssembler.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NavierStokesAssembler.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockAssembler.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/InterfaceAssembler.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NavierStokesStabilization.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SUPGStabilization.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/VMSStabilization.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/HFStabilization.hpp
    CACHE INTERNAL ""
)

set(SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/StokesAssembler_fem.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/StokesAssembler_rb.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockAssembler_fem.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockAssembler_rb.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/InterfaceAssembler_fem.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/InterfaceAssembler_rb.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NavierStokesAssembler_fem.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NavierStokesAssembler_rb.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/InletInflowAssembler_dense.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/InletInflowAssembler_sparse.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NavierStokesStabilization.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SUPGStabilization.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/VMSStabilization.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/HFStabilization.cpp
    CACHE INTERNAL ""
)
