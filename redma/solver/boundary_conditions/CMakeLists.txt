set(HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/BCManager.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/WindkesselModel.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PressureDrop.hpp
    CACHE INTERNAL ""
)

set(SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/BCManager.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/WindkesselModel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PressureDrop.cpp
    CACHE INTERNAL ""
)
