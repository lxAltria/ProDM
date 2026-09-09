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

${BUILD}/test_qoi_refactor data refactor 60 3 100 500 500 -d 0 0 > qoi.log
${BUILD}/test_qoi_reconstructor data refactor 1 0.01 -d 0 0 1 >> qoi.log

${BUILD}/test_pdr_refactor data/VelocityX.dat refactored 30 3 100 500 500 -d 4 > pdr.log
${BUILD}/test_pdr_reconstructor data/VelocityX.dat refactored 3 0.01 0.001 0.0001 -d 4 >> pdr.log

${BUILD}/test_qoi_refactor data refactor 60 3 100 500 500 -d 1 0 0.001 7 4 > qoi_weighted.log
${BUILD}/test_qoi_reconstructor data refactor 1 0.01 -d 1 0 1 >> qoi_weighted.log

${BUILD}/test_proaicd_refactor data/VelocityX.dat refactored -d 4 60 3 100 500 500 -Nega -eb -CP > proaicd.log
${BUILD}/test_proaicd_reconstructor data/VelocityX.dat refactored -d 3 0.01 0.001 0.0001 -Nega -DP -CP >> proaicd.log
