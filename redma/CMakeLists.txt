add_subdirectory(geometry)
add_subdirectory(solver)
add_subdirectory(utils)
add_subdirectory(reduced_basis)

# note: setting this to shared may cause problem at runtime when loading meshes
add_library(core STATIC
            ${HEADERS}
            ${SOURCES})

target_link_libraries(core ${LINKED_LIBRARIES})
