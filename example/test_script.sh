#!/bin/bash
# Run from the example folder: sh test_script.sh
# Requires the example data under data/ (see download_data.sh) and binaries under ../build
# Float data is kept as data/*.dat.f32; double data (data/*.dat) is derived by float2double.py
# The results are written to *.log

BUILD=../build/test
R=${BUILD}/prodm_refactor
Q=${BUILD}/prodm_retrieve

# multilevel decomposition [SC'21]
$R data refactored --vars VelocityX --dims 100 500 500 --dtype f --method mdr --target-level 4 --bitplanes 30 > mdr.log
$Q data refactored --vars VelocityX --dtype f --method mdr --tolerance 0.01 0.001 0.0001 >> mdr.log

# residual snapshots on SZ3 [TVCG'23]
$R data refactored --vars VelocityX --dims 100 500 500 --dtype f --method pdr-delta --approximator sz3 > pdr_delta.log
$Q data refactored --vars VelocityX --dtype f --method pdr-delta --approximator sz3 --tolerance 0.05 0.005 0.0005 >> pdr_delta.log

# convert the velocities to double (data/*.dat) for the QoI and PDR tests
if ! python3 -c "import numpy" 2>/dev/null; then
    echo "Error: python3 with numpy is required for float2double.py; activate an environment that provides it and re-run"
    exit 1
fi
for var in VelocityX VelocityY VelocityZ; do
    python3 float2double.py data/${var}.dat.f32
done

# QoI error control on V_total with the hand-derived estimator [SC'24]
$R data refactor --vars VelocityX,VelocityY,VelocityZ --dims 100 500 500 --dtype d --method pdr --approximator hpez > qoi.log
$Q data refactor --dtype d --method pdr --approximator hpez --qoi Vtot --estimator hand --inverse coordinate --joint-range --tolerance 0.01 >> qoi.log

# approximation-based refactoring with HPEZ [HPDC'26]
$R data refactored --vars VelocityX --dims 100 500 500 --dtype d --method pdr --approximator hpez --bitplanes 30 > pdr.log
$Q data refactored --vars VelocityX --dtype d --method pdr --approximator hpez --tolerance 0.01 0.001 0.0001 >> pdr.log

# QoI-weighted bitplanes [HPDC'26]
$R data refactor --vars VelocityX,VelocityY,VelocityZ --dims 100 500 500 --dtype d --method pdr --approximator hpez --weights hand:Vtot --eb 0.001 --max-weight 7 --block-size 4 > qoi_weighted.log
$Q data refactor --dtype d --method pdr --approximator hpez --qoi Vtot --estimator hand --inverse coordinate --joint-range --tolerance 0.01 >> qoi_weighted.log

# adaptive interpolation and coefficient decomposition [SC'26]
$R data refactored --vars VelocityX --dims 100 500 500 --dtype d --method proaicd --target-level 4 --bitplanes 60 --encoder nega --cp > proaicd.log
$Q data refactored --vars VelocityX --dtype d --method proaicd --encoder nega --interpreter dp --cp --tolerance 0.01 0.001 0.0001 >> proaicd.log
