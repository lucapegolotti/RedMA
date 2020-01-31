set(HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/AssemblerFactory.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/aAssembler.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/StokesAssembler.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NavierStokesAssembler.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockAssembler.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/InterfaceAssembler.hpp
    CACHE INTERNAL ""
)

set(SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/StokesAssembler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NavierStokesAssembler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/InterfaceAssembler.cpp
    CACHE INTERNAL ""
)
