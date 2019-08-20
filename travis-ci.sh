#!/bin/bash

set -e
set -u
set -x

sourcedir="$(pwd)"

cd ..

builddir="$(pwd)/cvc_depgraph_builddir"

###############################################################################
#                              Install Packages                               #
###############################################################################

ubuntu_packages=(
    cmake
    libhdf5-dev 
    hdf5-tools
    libeigen3-dev
    libboost-filesystem-dev libboost-system-dev libboost-program-options-dev libboost-graph-dev
    libgtest-dev
    openmpi-bin openmpi-common libopenmpi-dev
    libopenblas-base libopenblas-dev
    numdiff
)
sudo add-apt-repository "deb http://archive.ubuntu.com/ubuntu $(lsb_release -sc) main universe restricted multiverse"
sudo apt-get update
sudo apt-get install -y "${ubuntu_packages[@]}"

# Installation path for the various libraries to be installed below
mkdir -p local
install_prefix="$(realpath local)"

###############################################################################
#                               Install C-LIME                                #
###############################################################################
git clone https://github.com/usqcd-software/c-lime.git
pushd c-lime
./autogen.sh
./configure --prefix="$install_prefix"
make -j $(nproc)
make install
popd

export LIMEDIR="$install_prefix"

###############################################################################
#                               Install tmLQCD                                #
###############################################################################
git clone https://github.com/etmc/tmlqcd.git tmlqcd -b quda_work
pushd tmlqcd
tmlqcd_srcdir="$(pwd)"
autoconf
popd

rm -f tmlqcd_builddir
mkdir -p tmlqcd_builddir
pushd tmlqcd_builddir
CC=$(which mpicc) \
CFLAGS="-O3 -std=c99" \
"$tmlqcd_srcdir"/configure --disable-omp --enable-mpi --with-mpidimension=4 \
  --disable-sse2 --disable-sse3 \
  --enable-halfspinor --enable-gaugecopy \
  --enable-alignment=32 \
  --with-limedir="$install_prefix" \
  --with-lapack="-lblas -llapack"
make -j $(nproc)
tmlqcd_builddir="$(pwd)"
popd

###############################################################################
#                         Fetch and install HighFive lib                      #
###############################################################################
# note we use v2.0 as 'known-good'
highfive_builddir=highfive_builddir
git clone https://github.com/BlueBrain/HighFive.git -b v2.0
mkdir -p "$highfive_builddir"
pushd "$highfive_builddir"
cmake \
  -DCMAKE_INSTALL_PREFIX="$install_prefix" \
  -DHIGHFIVE_UNIT_TESTS=FALSE \
  -DHIGHFIVE_EXAMPLES=FALSE \
  ../HighFive 
make -j $(nproc)
make install
popd

###############################################################################
#                         Fetch and install yaml-cpp                          #
###############################################################################
yamlcpp_builddir=yamlcpp_builddir
git clone https://github.com/jbeder/yaml-cpp
mkdir -p "$yamlcpp_builddir"
pushd "$yamlcpp_builddir"
cmake \
  -DCMAKE_INSTALL_PREFIX="$install_prefix" \
  -DYAML_CPP_BUILD_TESTS=OFF \
  ../yaml-cpp
make -j $(nproc)
make install
popd

###############################################################################
#                         Build correlators executable                        #
###############################################################################
rm -rf "$builddir"
mkdir -p "$builddir"
pushd "$builddir"
CXX=$(which mpicxx)

cmake \
  -DCMAKE_CXX_FLAGS="-O3" \
  -DCMAKE_PREFIX_PATH="$install_prefix" \
  -DCMAKE_CXX_COMPILER="$CXX" \
  -DLIME_HOME="$install_prefix" \
  -DTMLQCD_SRC="$tmlqcd_srcdir" \
  -DTMLQCD_BUILD="$tmlqcd_builddir" \
  -DPARALLEL_LEVEL=TXYZ \
  "$sourcedir"

make -j $(nproc) || VERBOSE=1 make correlators

# set up environment variable that can be used in the integration test(s)
# to refer to correlators executable
export CORRBIN="$builddir/correlators"
popd

###############################################################################
#                         run integration test                                #
###############################################################################
pushd "$sourcedir"/integration_tests/correlators_8c16
./run_integration_test.sh
popd

##ctest --output-on-failure
