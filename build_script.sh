#!/bin/bash

source_dir=`pwd`
external_dir=${source_dir}/external
mkdir -p external
cd ${external_dir}
# build SZ (to use ZSTD compressor)
git clone https://github.com/szcompressor/SZ.git
cd SZ
git reset --hard f48d2f27a5470a28e900db9b46bb3344a2bc211f
mkdir -p build
mkdir -p install
cd build
cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_INSTALL_PREFIX=${external_dir}/SZ/install ..
make -j 8
make install

# build SZ3 (to use quantizer and huffman encoder)
cd ${external_dir}
git clone https://github.com/szcompressor/SZ3.git
cd SZ3
git reset --hard 89df80e547d839049202deb9d2bb9a11f4404bbb
mkdir -p build
cd build
cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_INSTALL_PREFIX=${external_dir}/SZ3/install ..
make -j 8
make install

# build MGARDx
cd ${external_dir}
git clone https://github.com/Linus-Li-1037/MGARDx.git
cd MGARDx
mkdir -p build
mkdir -p install
cd build
cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_INSTALL_PREFIX=${external_dir}/MGARDx/install ..
make -j 8
make install

# build ADIOS2
cd ${external_dir}
git clone https://github.com/ornladios/ADIOS2.git
cd ADIOS2
mkdir -p adios2-build && mkdir -p adios2-install
cd adios2-build
cmake -DADIOS2_USE_MPI=OFF -DADIOS2_USE_SZ=OFF -DCMAKE_INSTALL_PREFIX=${external_dir}/ADIOS2/adios2-install ..
make -j 8
make install

# build ProDM
cd ${source_dir}
mkdir -p build
cd build
cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ ..
make -j 8
