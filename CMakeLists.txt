CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

set (CMAKE_CXX_STANDARD 11)

project(RedMA VERSION 1.0.0 DESCRIPTION "Reduced Modeling of Arteries")

set(LINKED_LIBRARIES )

set(HEADERS_DIRS
    ${CMAKE_SOURCE_DIR}/geometric_structure
    ${CMAKE_SOURCE_DIR}/solver
    ${CMAKE_SOURCE_DIR}/utils
    ${TPL_MPI_INCLUDE_DIRS}
    ${TPL_Trilinos_INCLUDE_DIRS}
    ${TPL_LifeV_INCLUDE_DIRS}
    ${TPL_HDF5_INCLUDE_DIRS}
    ${TPL_BOOST_INCLUDE_DIRS}
    ${TPL_PARMETIS_INCLUDE_DIRECTORY}
    ${TPL_METIS_INCLUDE_DIRECTORY}
    ${TPL_SUITESPARSE_INCLUDE_DIRECTORY}
    ${TPL_OPENBLAS_INCLUDE_DIRECTORY}
    ${TPL_TINYXML2_INCLUDE_DIRECTORY}
    )

set(LIBRARIES_DIRS
    ${TPL_MPI_LIBRARY_DIRS}
    ${TPL_Trilinos_LIBRARY_DIRS}
    ${TPL_LifeV_LIBRARY_DIRS}
    ${TPL_HDF5_LIBRARY_DIRS}
    ${TPL_BOOST_LIBRARY_DIRS}
    ${TPL_PARMETIS_LIBRARY_DIRECTORY}
    ${TPL_METIS_LIBRARY_DIRECTORY}
    ${TPL_SUITESPARSE_LIBRARY_DIRECTORY}
    ${TPL_OPENBLAS_LIBRARY_DIRECTORY}
    ${TPL_TINYXML2_LIBRARY_DIRECTORY}
    )

set(TRILINOS_LIBRARIES
    amesos
    anasazi
    anasaziepetra
    aztecoo
    belos
    belosepetra
    epetra
    epetraext
    galeri-epetra
    ifpack
    isorropia
    kokkoscore
    ml
    ModeLaplace
    rtop
    teuchoscomm
    teuchoscore
    teuchoskokkoscomm
    teuchoskokkoscompat
    teuchosnumerics
    teuchosparameterlist
    teuchosremainder
    thyracore
    thyraepetra
    thyraepetraext
    trilinosss
    triutils
    zoltan
    belosepetra)

set(LINKED_LIBRARIES
    boost_system
    boost_filesystem
    tinyxml2
    mpi
    mpi_cxx
    lifevcore
    lifeveta
    lifevnavierstokesblocks
    ${TRILINOS_LIBRARIES}
    hdf5
    parmetis
    metis
    lapack
    umfpack
    openblas
    amd
    )

include_directories(${HEADERS_DIRS})
link_directories(${LIBRARIES_DIRS})

add_subdirectory(geometric_structure)
add_subdirectory(solver)
add_subdirectory(utils)
add_subdirectory(examples)
add_subdirectory(tests)
