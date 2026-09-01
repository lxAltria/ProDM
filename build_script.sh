#!/bin/bash
set -eu

export CMAKE_POLICY_VERSION_MINIMUM=3.5
# Compilers: uses the system default unless overridden, e.g.
#   CC=gcc-16 CXX=g++-16 sh build_script.sh
# Note: changing compilers requires removing the build directories first,
# since CMake caches the compiler on the first configure.
export CC=${CC:-gcc}
export CXX=${CXX:-g++}

# Pinned dependency versions (tested together; QoZ headers must match SZ3's)
SZ2_COMMIT=308bd06f0040ec0d5c22fb3fcb0428c306ba4df1
SZ3_COMMIT=90c66bed1c04e701442ecb104b912548fcfabee9
QOZ_COMMIT=17f124ef9a341f4bd661c20e6f0a3e80b19c43ef
MGARD_COMMIT=7ba6738d429a70da8bd1d345ac7aad702366f16b

# clone (if needed) and check out the pinned commit
fetch_dependency() {
    name=$1; url=$2; commit=$3
    if [ ! -d "${name}" ]; then
        git clone "${url}" "${name}"
    fi
    cd "${name}"
    git cat-file -e "${commit}^{commit}" 2>/dev/null || git fetch origin
    git reset --hard "${commit}"
    cd ..
}

source_dir=`pwd`
external_dir=${source_dir}/external
mkdir -p external
cd ${external_dir}

# build SZ2
fetch_dependency SZ2 https://github.com/szcompressor/SZ2.git ${SZ2_COMMIT}
cd SZ2
mkdir -p build
mkdir -p install
cd build
cmake -DCMAKE_INSTALL_PREFIX=${external_dir}/SZ2/install ..
make -j 8
make install

# build SZ3
cd ${external_dir}
fetch_dependency SZ3 https://github.com/szcompressor/SZ3.git ${SZ3_COMMIT}
cd SZ3
mkdir -p build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${external_dir}/SZ3/install ..
make -j 8
make install

# build QoZ (HPEZ) -- header-only; ska_hash is excluded because QoZ reuses SZ3's copy
cd ${external_dir}
fetch_dependency QoZ https://github.com/Linus-Li-1037/QoZ.git ${QOZ_COMMIT}
cd QoZ
mkdir -p install/include/
rsync -a --exclude 'ska_hash' include/ ${external_dir}/QoZ/install/include/

# build MGARD
cd ${external_dir}
fetch_dependency MGARD https://github.com/CODARcode/MGARD.git ${MGARD_COMMIT}
cd MGARD
sh build_scripts/build_mgard_serial.sh 8

# build ProDM
cd ${source_dir}
mkdir -p build
mkdir -p install
cd build
cmake -DCMAKE_INSTALL_PREFIX=${source_dir}/install -DPRODM_WITH_SZ2=ON -DPRODM_WITH_SZ3=ON -DPRODM_WITH_HPEZ=ON -DPRODM_WITH_MGARD=ON ..
make -j 8

# build ADIOS2
# cd ${external_dir}
# git clone https://github.com/ornladios/ADIOS2.git
# cd ADIOS2
# mkdir -p adios2-build && mkdir -p adios2-install
# cd adios2-build
# cmake -DADIOS2_USE_MPI=OFF -DADIOS2_USE_SZ=OFF -DCMAKE_INSTALL_PREFIX=${external_dir}/ADIOS2/adios2-install ..
# make -j 8
# make install
