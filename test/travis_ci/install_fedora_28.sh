#!/bin/bash
set -v

# suck in package lists
dnf update -qq -y

# install deps
# use PIP for Python packages
dnf install -qq -y environment-modules which git-all gcc-c++ gcc-gfortran \
    make cmake-3.11.0-1.fc28 swig mpich-devel hdf5-mpich-devel \
    netcdf-mpich-devel boost-devel python3-devel python3-pip subversion \
    udunits2 udunits2-devel zlib-devel openssl-devel wget

source /usr/share/Modules/init/bash

module load mpi/mpich-x86_64

echo ${TRAVIS_BRANCH}
echo ${BUILD_TYPE}
echo ${DOCKER_IMAGE}
echo ${IMAGE_VERSION}
echo ${TECA_PYTHON_VERSION}
echo ${TECA_DATA_REVISION}

pip${TECA_PYTHON_VERSION} install numpy mpi4py matplotlib torch

# install data files.
svn co svn://svn.code.sf.net/p/teca/TECA_data@${TECA_DATA_REVISION} TECA_data
