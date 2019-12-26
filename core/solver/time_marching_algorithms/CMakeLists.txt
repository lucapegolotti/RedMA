set(HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/BackwardEuler.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RosenbrockAlgorithm.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SteadySolver.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TimeMarchingAlgorithm.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TimeMarchingAlgorithmsFactory.hpp
    CACHE INTERNAL ""
)

set(SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/BackwardEuler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RosenbrockAlgorithm.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SteadySolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TimeMarchingAlgorithm.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TimeMarchingAlgorithmsFactory.cpp
    CACHE INTERNAL ""
)
