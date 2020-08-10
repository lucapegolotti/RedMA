set(HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/TimeMarchingAlgorithmFactory.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/aTimeMarchingAlgorithm.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BDF.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GeneralizedAlphaMethod.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/aFunctionProvider.hpp
    CACHE INTERNAL ""
)

set(SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/TimeMarchingAlgorithmFactory.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/aTimeMarchingAlgorithm.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BDF.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GeneralizedAlphaMethod.cpp
    CACHE INTERNAL ""
)
