module restore juwels_intelgccparastation \
module load CMake \
CXX=CC \
CC=cc \
CXXFLAGS="-fopenmp -O3 -mtune=haswell -march=haswell -g" \
CFLAGS="-fopenmp -O3 -mtune=haswell -march=haswell -g" \
LDFLAGS="-fopenmp" \
cmake \
  -DLIME_HOME=/p/home/jusers/pittler1/juwels/build/lime \
  -DLEMON_HOME=/p/project/chbn28/hbn28d/build/lemon \
/p/project/chbn28/hbn28d/code/cvc-1
