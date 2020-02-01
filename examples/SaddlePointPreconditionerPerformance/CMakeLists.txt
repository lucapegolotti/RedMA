add_executable(SaddlePointPreconditionerPerformance SaddlePointPreconditionerPerformance.cpp)
target_link_libraries(SaddlePointPreconditionerPerformance ${LINKED_LIBRARIES})
target_include_directories(SaddlePointPreconditionerPerformance PRIVATE ${HEADERS_DIRS})

add_custom_command(
    TARGET SaddlePointPreconditionerPerformance PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/datafiles
    ${CMAKE_CURRENT_BINARY_DIR}/datafiles
)

add_custom_command(
    TARGET SaddlePointPreconditionerPerformance PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/geometries
    ${CMAKE_CURRENT_BINARY_DIR}/geometries
)

add_custom_command(
    TARGET SaddlePointPreconditionerPerformance PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/scripts
    ${CMAKE_CURRENT_BINARY_DIR}/scripts
)
