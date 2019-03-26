#!/bin/bash
. /users/fpittler/build/cvc2/module_loads.sh
CXX=CC \
CC=cc \
CXXFLAGS="-fopenmp -O3 -mtune=haswell -march=haswell -g" \
CFLAGS="-fopenmp -O3 -mtune=haswell -march=haswell -g" \
LDFLAGS="-fopenmp" \
cmake \
  -DLAPACK_HOME=/opt/cray/pe/libsci/18.07.1/GNU/6.1/x86_64 \
  -DLAPACK_LIBRARIES="-lsci_gnu_mp" \
  -DBLAS_HOME=/opt/cray/pe/libsci/18.07.1/GNU/6.1/x86_64 \
  -DBLAS_LIBRARIES="-lsci_gnu_mp" \
  -DLIME_HOME=/users/bartek/local/haswell/libs/lime \
  -DLEMON_HOME=/users/bartek/local/haswell/libs/lemon \
 /users/fpittler/code/cvc_new/cvc-1
