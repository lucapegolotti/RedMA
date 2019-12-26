set(HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockMatrix.hpp
    CACHE INTERNAL ""
)

set(SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockMatrix.cpp
    CACHE INTERNAL ""
)

add_library(array SHARED
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockMatrix.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BlockMatrix.cpp
)

target_link_libraries(array ${LINKED_LIBRARIES} ${HEADERS} ${SOURCES})
