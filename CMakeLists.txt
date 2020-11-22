CMAKE_MINIMUM_REQUIRED(VERSION 3.14)

project(RedMA VERSION 1.0.0 DESCRIPTION "Reduced Modeling of Arteries")

set(CMAKE_CXX_STANDARD 17)

set(SOURCES "" CACHE INTERNAL "")
set(HEADERS "" CACHE INTERNAL "")

# message(STATUS ${TPL_LifeV_INCLUDE_DIRS})
set(HEADERS_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
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
    ${TPL_PYTHON_INCLUDE_DIRECTORY}
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
    belosepetra
    ${TPL_PYTHON_LIBRARY}
    )

# note: by default, boost is built without additional libraries in the lifev
# environment. If that version of boost is used, boost must be installed by removing
# the -without-system and -without-filesystem flags in the bootstrap command of lifev
set(LINKED_LIBRARIES
    tinyxml2
    mpi
    lifevcore
    lifeveta
    lifevnavierstokesblocks
    #  lifevreducedbasis
    ${TRILINOS_LIBRARIES}
    hdf5
    parmetis
    metis
    umfpack
    openblas
    amd
    )

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# we disable warnings of wrong return type because these are already present
# at the level of LifeV
set(CMAKE_CXX_FLAGS "-std=c++17 -Wno-return-type -Wno-unknown-pragmas")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

link_directories(${LIBRARIES_DIRS})
include_directories(${HEADERS_DIRS})

# this is to avoid problems of rpath when launching targets linking to core
# set(CMAKE_SKIP_BUILD_PATH FALSE)
# set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
# set(CMAKE_INSTALL_RPATH ${LIBRARIES_DIRS})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

add_subdirectory(redma)

# add library libcore to linked libraries
set(LINKED_LIBRARIES ${LINKED_LIBRARIES} core)

add_subdirectory(examples)
enable_testing()
add_subdirectory(tests)
