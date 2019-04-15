file(COPY
     ${CMAKE_CURRENT_SOURCE_DIR}/bifurcation_symmetric_coarse.mesh
     DESTINATION
     ${CMAKE_CURRENT_BINARY_DIR}
     )

file(COPY
     ${CMAKE_CURRENT_SOURCE_DIR}/tube_coarse.mesh
     DESTINATION
     ${CMAKE_CURRENT_BINARY_DIR}
     )

file(COPY
     ${CMAKE_CURRENT_SOURCE_DIR}/datafiles
     DESTINATION
     ${CMAKE_CURRENT_BINARY_DIR}
     )
