#!/bin/bash

export CMAKE_POLICY_VERSION_MINIMUM=3.5
export CC=gcc-16
export CXX=g++-16
source_dir=`pwd`
external_dir=${source_dir}/external
mkdir -p external
cd ${external_dir}

# build SZ2
git clone https://github.com/szcompressor/SZ2.git
cd SZ2
mkdir -p build
mkdir -p install
cd build
cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_INSTALL_PREFIX=${external_dir}/SZ2/install ..
make -j
make install

# build SZ3
cd ${external_dir}
git clone https://github.com/szcompressor/SZ3.git
cd SZ3
git reset --hard 90c66bed1c04e701442ecb104b912548fcfabee9
mkdir -p build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${external_dir}/SZ3/install ..
make -j
make install

# build MGARDx
cd ${external_dir}
git clone https://github.com/Linus-Li-1037/MGARDx.git
cd MGARDx
mkdir -p build
mkdir -p install
cd build
cmake -DCMAKE_INSTALL_PREFIX=${external_dir}/MGARDx/install ..
make -j
make install

# build QoZ (HPEZ)
cd ${external_dir}
git clone https://github.com/Linus-Li-1037/QoZ.git
cd QoZ
mkdir -p install/include/
rsync -av --exclude 'ska_hash' include/ ${external_dir}/QoZ/install/include/

# build MGARD
cd ${external_dir}
git clone https://github.com/CODARcode/MGARD.git
cd MGARD
sh build_scripts/build_mgard_serial.sh 8

# build ProDM
cd ${source_dir}
mkdir -p build
mkdir -p install
cd build
cmake -DCMAKE_INSTALL_PREFIX=${source_dir}/install ..
make -j 8

# build ADIOS2
cd ${external_dir}
git clone https://github.com/ornladios/ADIOS2.git
cd ADIOS2
mkdir -p adios2-build && mkdir -p adios2-install
cd adios2-build
cmake -DADIOS2_USE_MPI=OFF -DADIOS2_USE_SZ=OFF -DCMAKE_INSTALL_PREFIX=${external_dir}/ADIOS2/adios2-install ..
make -j 8
make install

