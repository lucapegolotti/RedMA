add_subdirectory(abstract)
add_subdirectory(block)
add_subdirectory(coupling)
add_subdirectory(finite_element)
add_subdirectory(reduced_basis)
add_subdirectory(models)

file(GLOB headers_files "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
set(HEADERS ${HEADERS} ${headers_files} CACHE INTERNAL "")

file(GLOB source_files "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
set(SOURCES ${SOURCES} ${source_files} CACHE INTERNAL "")
