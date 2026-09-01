#!/bin/bash

mkdir -p data
curl -LO https://g-d0cd3f.fd635.8443.data.globus.org/raw-data/Hurricane-ISABEL/SDRBENCH-Hurricane-ISABEL-100x500x500.tar.gz
tar -xvf SDRBENCH-Hurricane-ISABEL-100x500x500.tar.gz
cp 100x500x500/Uf48.bin.f32 data/VelocityX.dat
cp 100x500x500/Vf48.bin.f32 data/VelocityY.dat
cp 100x500x500/Wf48.bin.f32 data/VelocityZ.dat
