file(GLOB headers_files "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
set(HEADERS ${HEADERS} ${headers_files} CACHE INTERNAL "")

file(GLOB source_files "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
set(SOURCES ${SOURCES} ${source_files} CACHE INTERNAL "")
