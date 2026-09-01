#!/bin/bash

git checkout ../example/data/VelocityX.dat
./test/test_mdr_refactor ../example/data/VelocityX.dat refactored 4 30 3 100 500 500 0 -f > mdr.log
./test/test_mdr_reconstructor ../example/data/VelocityX.dat refactored 3 0.01 0.001 0.0001 0 -f >> mdr.log

./test/test_pdr_delta_refactor ../example/data/VelocityX.dat refactored 3 100 500 500 -f 3 > pdr_delta.log
./test/test_pdr_delta_reconstructor ../example/data/VelocityX.dat refactored 3 0.05 0.005 0.0005 -f 3 >> pdr_delta.log

source /Users/xin/py_env/bin/activate
python3 ../example/float2double.py ../example/data/VelocityX.dat
deactivate

./test/refactor_d64 4 0 Hurricane ../example > qoi.log
./test/qoi_Vtot_d64 4 0 1 0.01 ../example >> qoi.log

./test/test_pdr_refactor ../example/data/VelocityX.dat refactored 30 3 100 500 500 -d 4 > pdr.log
./test/test_pdr_reconstructor ../example/data/VelocityX.dat refactored 3 0.01 0.001 0.0001 -d 4 >> pdr.log

./test/refactor_d64 4 1 Hurricane ../example 7 4 0.001 > qoi_weighted.log
./test/qoi_Vtot_d64 4 1 1 0.01 ../example >> qoi_weighted.log

./test/two_modes_refactor ../example/data/VelocityX.dat refactored -d 4 60 3 100 500 500 -Nega -eb -CP > proaicd.log
./test/two_modes_reconstructor ../example/data/VelocityX.dat refactored -d 3 0.01 0.001 0.0001 -Nega -DP -CP >> proaicd.log

