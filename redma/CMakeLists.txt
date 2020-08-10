add_subdirectory(array)
add_subdirectory(assemblers)
add_subdirectory(boundary_conditions)
add_subdirectory(coupling_basis_functions)
add_subdirectory(problem)
add_subdirectory(reduced_basis)
add_subdirectory(solver)
add_subdirectory(utils)

# note: setting this to shared may cause problem at runtime when loading meshes
add_library(core STATIC
            ${HEADERS}
            ${SOURCES})

target_link_libraries(core ${LINKED_LIBRARIES})
