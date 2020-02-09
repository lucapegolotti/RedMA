set(HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/TimeMarchingAlgorithmFactory.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TimeMarchingAlgorithmFactory_imp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/aTimeMarchingAlgorithm.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/aTimeMarchingAlgorithm_imp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BDF.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BDF_imp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GeneralizedAlphaMethod.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GeneralizedAlphaMethod_imp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GeneralizedAlphaMethod1stOrderPressure.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GeneralizedAlphaMethod1stOrderPressure_imp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/aFunctionProvider.hpp
    CACHE INTERNAL ""
)

set(SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/GeneralizedAlphaMethod1stOrderPressure.cpp
    CACHE INTERNAL ""
)
