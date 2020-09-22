#!/bin/bash

#export CC=gcc
#export CXX=g++
#export CMAKE_C_COMPILER=$CC
#export CMAKE_CXX_COMPILER=$CXX
#export CMAKE_C_LINK_EXECUTABLE=$CC
#export CMAKE_CXX_LINK_EXECUTABLE=$CXX

export SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export COVERAGE="-g -O0 -fprofile-arcs -ftest-coverage"
export COVERAGE=""
export BUILD_DIR="$SRC_DIR/build"
export INSTALL_DIR=$HOME/apbs
export PATH=$INSTALL_DIR:$PATH
export RELEASE_TYPE=Debug
export RELEASE_TYPE=Release

echo "==================================== WHERE AM I ==================================== "
pwd
echo "==================================== VERSIONS: ==================================== "
echo "==================================== PYTHON VERSION"
python -c "import sys; print(sys.version)"
echo "==================================== CMAKE VERSION"
cmake --version
echo "==================================== C Compiler VERSION"
$CMAKE_C_COMPILER --version
echo "==================================== C++ Compiler VERSION"
$CMAKE_CXX_COMPILER --version
echo "==================================== SWIG VERSION"
swig -version
echo "==================================== Install Python requirements ==================================== "
pip3 install -U pip
pip3 install -U pytest
pip3 install -U virtualenv
pip3 install -U numpy
pip3 install -r requirements.txt
#  Just build APBS for now
echo "==================================== PWD FOR TOP DIR ==================================== "
pwd
echo "==================================== Get External SubModules ==================================== "
git submodule init
git submodule update
echo "==================================== BUILD =============================================== "

rm -rf $BUILD_DIR                                         || exit 1
rm -rf $INSTALL_DIR                                       || exit 1
mkdir -p $BUILD_DIR                                       || exit 1
mkdir -p $INSTALL_DIR                                     || exit 1
#  Build pybind11
pushd $(pwd)/externals/pybind11
[ -d build ] || mkdir -p build
[ -d install ] || mkdir -p install
pushd build
cmake .. -DDOWNLOAD_CATCH=ON -DCMAKE_INSTALL_PREFIX=$(python3 -c "import os,sys; print(os.path.realpath(sys.argv[1]))" ../install)
make -j install
popd
export pybind11_DIR=$(python3 -c "import os,sys; print(os.path.realpath(sys.argv[1]))" ./install)
popd
cd $BUILD_DIR                                             || exit 1
#cmake -S .. -B $BUILD_DIR --trace-source=../CMakeLists.txt --trace-expand \
cmake                                                     \
      -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR                 \
      -DCMAKE_BUILD_TYPE=$RELEASE_TYPE                    \
      -DENABLE_GEOFLOW=ON                                 \
      -DENABLE_BEM=ON                                     \
      -DENABLE_FETK=ON                                    \
      -DENABLE_OPENMP=ON                                  \
      -DENABLE_PBAM=ON                                    \
      -DENABLE_PBSAM=ON                                   \
      -DENABLE_PYTHON=ON                                  \
      -DENABLE_TESTS=ON                                   \
      -DENABLE_TINKER=OFF                                 \
      -DBUILD_SHARED_LIBS=ON                              \
      -DBUILD_DOC=OFF                                     \
      -DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} ${COVERAGE}"      \
      -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} ${COVERAGE}"  \
      ..                                                  || exit 1
#      -DCMAKE_C_FLAGS="-fPIC"                             \
VERBOSE=1 make -j 1                                       || exit 1
VERBOSE=1 make -j 1 install                               #|| exit 1
export PATH="$INSTALL_DIR/bin:$PATH"
# Run a single test if it fails using the following:
# ctest -VV -R pbam_test
ctest -C Release --output-on-failure                      #|| exit 1

cpack -C Release -G ZIP                                   || exit 1
unzip -l APBS*.zip
