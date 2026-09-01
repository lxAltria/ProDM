#!/bin/bash
# Run from the example folder: sh test_script.sh
# Requires the example data under data/ (see download_data.sh) and binaries under ../build
# Float data is kept as data/*.dat.f32; double data (data/*.dat) is derived by float2double.py

BUILD=../build/test

${BUILD}/test_mdr_refactor data/VelocityX.dat.f32 refactored 4 30 3 100 500 500 0 -f > mdr.log
${BUILD}/test_mdr_reconstructor data/VelocityX.dat.f32 refactored 3 0.01 0.001 0.0001 0 -f >> mdr.log

${BUILD}/test_pdr_delta_refactor data/VelocityX.dat.f32 refactored 3 100 500 500 -f 3 > pdr_delta.log
${BUILD}/test_pdr_delta_reconstructor data/VelocityX.dat.f32 refactored 3 0.05 0.005 0.0005 -f 3 >> pdr_delta.log

# convert the velocities to double (data/*.dat) for the QoI and PDR tests
if ! python3 -c "import numpy" 2>/dev/null; then
    echo "Error: python3 with numpy is required for float2double.py; activate an environment that provides it and re-run"
    exit 1
fi
for var in VelocityX VelocityY VelocityZ; do
    python3 float2double.py data/${var}.dat.f32
done

${BUILD}/refactor_d64 4 0 Hurricane . > qoi.log
${BUILD}/qoi_Vtot_d64 4 0 1 0.01 . >> qoi.log

${BUILD}/test_pdr_refactor data/VelocityX.dat refactored 30 3 100 500 500 -d 4 > pdr.log
${BUILD}/test_pdr_reconstructor data/VelocityX.dat refactored 3 0.01 0.001 0.0001 -d 4 >> pdr.log

${BUILD}/refactor_d64 4 1 Hurricane . 7 4 0.001 > qoi_weighted.log
${BUILD}/qoi_Vtot_d64 4 1 1 0.01 . >> qoi_weighted.log

${BUILD}/two_modes_refactor data/VelocityX.dat refactored -d 4 60 3 100 500 500 -Nega -eb -CP > proaicd.log
${BUILD}/two_modes_reconstructor data/VelocityX.dat refactored -d 3 0.01 0.001 0.0001 -Nega -DP -CP >> proaicd.log
