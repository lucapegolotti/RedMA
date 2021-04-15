#!/bin/bash
module load cmake
module load gcc/7.4.0
module load openblas
module load mvapich2

THIS=${PWD}

rm CMakeCache.txt
rm -rf CMakeFiles
mkdir /scratch/nicmuell/RedMA_build/build

# if you have installed LifeV by following the instructions on www.lifev.org,
# LIBRARIES_BASE_DIRECTORY is a sub-directory of your lifev-env folder

LIBRARIES_BASE_DIRECTORY=/home/nicmuell/libraries/lifev-env/libs/installs/

HDF5_INCLUDE_DIR=${LIBRARIES_BASE_DIRECTORY}hdf5-1.8.19_installRelease/include/
HDF5_LIB_DIR=${LIBRARIES_BASE_DIRECTORY}hdf5-1.8.19_installRelease/lib/
TRILINOS_INCLUDE_DIR=${LIBRARIES_BASE_DIRECTORY}trilinos-12.12.1_installRelease/include/
TRILINOS_LIB_DIR=${LIBRARIES_BASE_DIRECTORY}trilinos-12.12.1_installRelease/lib/
BOOST_LIB_DIR=${LIBRARIES_BASE_DIRECTORY}boost_1_65_1_installRelease/lib/
BOOST_INCLUDE_DIR=${LIBRARIES_BASE_DIRECTORY}boost_1_65_1_installRelease/include/
PARMETIS_INCLUDE_DIRECTORY=${LIBRARIES_BASE_DIRECTORY}parmetis-4.0.3_installRelease/include/
PARMETIS_LIB_DIRECTORY=${LIBRARIES_BASE_DIRECTORY}parmetis-4.0.3_installRelease/lib/
METIS_INCLUDE_DIRECTORY=${LIBRARIES_BASE_DIRECTORY}metis-5.1.0_installRelease/include/
METIS_LIB_DIRECTORY=${LIBRARIES_BASE_DIRECTORY}metis-5.1.0_installRelease/lib/
SUITESPARSE_INCLUDE_DIRECTORY=${LIBRARIES_BASE_DIRECTORY}suitesparse-4.5.1_installRelease/include/
SUITESPARSE_LIB_DIRECTORY=${LIBRARIES_BASE_DIRECTORY}suitesparse-4.5.1_installRelease/lib/
# OPENBLAS_INCLUDE_DIRECTORY=${LIBRARIES_BASE_DIRECTORY}OpenBLAS-0.2.20_installRelease/include/
OPENBLAS_INCLUDE_DIRECTORY=/ssoft/spack/humagne/v1/opt/spack/linux-rhel7-x86_S6g1_Mellanox/gcc-7.4.0/openblas-0.3.6-t3zk7snalnypxb3ee3i7rxrpvbacejto/include/
# OPENBLAS_LIB_DIRECTORY=${LIBRARIES_BASE_DIRECTORY}OpenBLAS-0.2.20_installRelease/lib/
# OPENBLAS_LIB_DIRECTORY=/ssoft/spack/humagne/v1/opt/spack/linux-rhel7-x86_S6g1_Mellanox/gcc-7.4.0/openblas-0.3.6-t3zk7snalnypxb3ee3i7rxrpvbacejto/lib/
OPENBLAS_LIB_DIRECTORY=/ssoft/spack/humagne/v1/opt/spack/linux-rhel7-x86_S6g1_Mellanox/gcc-7.4.0/openblas-0.3.6-t3zk7snalnypxb3ee3i7rxrpvbacejto/lib
TINYXML2_INCLUDE_DIRECTORY=/home/nicmuell/libraries/tinyxml2/
TINYXML2_LIB_DIRECTORY=/home/nicmuell/libraries/tinyxml2/
MPI_INCLUDE_DIRECTORY=/ssoft/spack/humagne/v1/opt/spack/linux-rhel7-x86_E5v4_Mellanox/gcc-7.4.0/mvapich2-2.3.1-nvnz73udfwp4zqg47ttdrpodvvc2zvc5/include/
MPI_LIB_DIRECTORY=/ssoft/spack/humagne/v1/opt/spack/linux-rhel7-x86_E5v4_Mellanox/gcc-7.4.0/mvapich2-2.3.1-nvnz73udfwp4zqg47ttdrpodvvc2zvc5/lib/
# MPI_INCLUDE_DIRECTORY=/ssoft/spack/humagne/v1/opt/spack/linux-rhel7-x86_S6g1_Mellanox/gcc-8.3.0/mvapich2-2.3.1-ka5klaihzoxvao3kxweyvxjkh2soydog/include
# MPI_LIB_DIRECTORY=/ssoft/spack/humagne/v1/opt/spack/linux-rhel7-x86_S6g1_Mellanox/gcc-8.3.0/mvapich2-2.3.1-ka5klaihzoxvao3kxweyvxjkh2soydog/lib

LIFEV_INSTALLATION=/home/nicmuell/libraries/lifev-env/lifev-epfl-install/

BUILD_FDR=/scratch/nicmuell/RedMA_build_fidis/
mkdir $BUILD_FDR 
cp -r meshes $BUILD_FDR
mkdir $BUILD_FDR/build
cd $BUILD_FDR/build

cmake \
-D CMAKE_BUILD_TYPE:STRING=RELEASE \
-D TPL_LifeV_INCLUDE_DIRS:PATH=${LIFEV_INSTALLATION}include/ \
-D TPL_LifeV_LIBRARY_DIRS:PATH=${LIFEV_INSTALLATION}lib/ \
-D TPL_HDF5_INCLUDE_DIRS:PATH=${HDF5_INCLUDE_DIR} \
-D TPL_HDF5_LIBRARY_DIRS:PATH=${HDF5_LIB_DIR} \
-D TPL_Trilinos_INCLUDE_DIRS:PATH=${TRILINOS_INCLUDE_DIR} \
-D TPL_Trilinos_LIBRARY_DIRS:PATH=${TRILINOS_LIB_DIR} \
-D TPL_BOOST_INCLUDE_DIRS:PATH=${BOOST_INCLUDE_DIR} \
-D TPL_BOOST_LIBRARY_DIRS:PATH=${BOOST_LIB_DIR} \
-D TPL_MPI_INCLUDE_DIRS:PATH=${MPI_INCLUDE_DIRECTORY} \
-D TPL_MPI_LIBRARY_DIRS:PATH=${MPI_LIB_DIRECTORY} \
-D TPL_PARMETIS_INCLUDE_DIRECTORY:PATH=${PARMETIS_INCLUDE_DIRECTORY} \
-D TPL_PARMETIS_LIBRARY_DIRECTORY:PATH=${PARMETIS_LIB_DIRECTORY} \
-D TPL_METIS_INCLUDE_DIRECTORY:PATH=${METIS_INCLUDE_DIRECTORY} \
-D TPL_METIS_LIBRARY_DIRECTORY:PATH=${METIS_LIB_DIRECTORY} \
-D TPL_SUITESPARSE_INCLUDE_DIRECTORY:PATH=${SUITESPARSE_INCLUDE_DIRECTORY} \
-D TPL_SUITESPARSE_LIBRARY_DIRECTORY:PATH=${SUITESPARSE_LIB_DIRECTORY} \
-D TPL_OPENBLAS_LIBRARY_DIRECTORY:PATH=${OPENBLAS_LIB_DIRECTORY} \
-D TPL_OPENBLAS_INCLUDE_DIRECTORY:PATH=${OPENBLAS_INCLUDE_DIRECTORY} \
-D TPL_TINYXML2_LIBRARY_DIRECTORY:PATH=${TINYXML2_LIB_DIRECTORY} \
-D TPL_TINYXML2_INCLUDE_DIRECTORY:PATH=${TINYXML2_INCLUDE_DIRECTORY} \
$THIS

make -j 8
