set(HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/BuildingBlock.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Tube.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Aorta.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BifurcationSymmetric.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TreeStructure.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GeometryParser.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NonAffineDeformer.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GeometricParametersHandler.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GeometryPrinter.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SegmentationParser.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SegmentationsMerger.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AortaBifurcation0.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AortaBifurcation1.hpp
    CACHE INTERNAL ""
)

set(SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/BuildingBlock.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Tube.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Aorta.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AortaBifurcation0.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AortaBifurcation1.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BifurcationSymmetric.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TreeStructure.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GeometryParser.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NonAffineDeformer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GeometricParametersHandler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GeometryPrinter.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SegmentationParser.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SegmentationsMerger.cpp
    CACHE INTERNAL ""
)
