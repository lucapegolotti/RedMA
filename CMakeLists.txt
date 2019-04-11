CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

project(ReMA VERSION 1.0.0 DESCRIPTION "Reduced Modeling of Arteries")

set(LINKED_LIBRARIES tinyxml2)

set(HEADERS_DIRS ${CMAKE_SOURCE_DIR}/geometric_structure
                 ${CMAKE_SOURCE_DIR}/utils)

add_subdirectory(geometric_structure)
add_subdirectory(utils)
add_subdirectory(examples)
