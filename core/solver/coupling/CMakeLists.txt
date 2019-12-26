set(HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/Coupler.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Interpolation.hpp
    CACHE INTERNAL ""
)

set(SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/Coupler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Interpolation.cpp
    CACHE INTERNAL ""
)

# add_library(coupling SHARED
#             ${CMAKE_CURRENT_SOURCE_DIR}/Coupler.cpp
#             ${CMAKE_CURRENT_SOURCE_DIR}/Interpolation.cpp)
