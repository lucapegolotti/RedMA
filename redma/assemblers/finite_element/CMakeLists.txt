set(HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/StokesAssemblerFE.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/HFStabilization.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NavierStokesStabilization.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SUPGStabilization.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/VMSStabilization.hpp
    CACHE INTERNAL ""
)

set(SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/StokesAssemblerFE.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/HFStabilization.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NavierStokesStabilization.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SUPGStabilization.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/VMSStabilization.cpp
    CACHE INTERNAL ""
)
