# ProDM: Progressive Data Refactoring and Retrieval Framework

This is the code repo for NSF project "Collaborative Research: Elements: ProDM: Developing A Unified Progressive Data Management Library for Exascale Computational Science". It is a joint collaborative effort from the Oregon State University (OSU), New Jersey Institute of Technology (NJIT), and Temple University.

## Authors

**Major contributors:**
Wenbo Li (OSU), Xuan Wu (OSU), Qirui Tian (NJIT)

**Supervisors:**
Dr. Xin Liang (OSU), Dr. Qing Liu (NJIT), Dr. Xubin He (Temple)

**Other contributors and collaborators:**
Dr. Scott Klasky (ORNL), Dr. Qian Gong (ORNL), Dr. Jieyang Chen (Univ. of Oregon), Dr. Jill Zhang (LLNL), Dr. Seung-Hoe Ku (PPPL), Dr. Xiaohua Zhang (LLNL), etc.

## Publications
ProDM hosts multiple novel progressive compression methods developed by the research group:
- **[SC'26]**: Wenbo Li, Xuan Wu, Qian Gong, Pu Jiao, Jieyang Chen, Qing Liu, Norbert Podhorszki, Scott Klasky, and Xin Liang. *Improving Progressive Compression with Adaptive Interpolation and Coefficient Decomposition*
- **[HPDC'26]**: Wenbo Li, Qian Gong, Xuan Wu, Jieyang Chen, Qing Liu, Xubin He, Norbert Podhorszki, Scott Klasky, and Xin Liang. *QProR: An Efficient Framework for Quantity-of-Interest Based Progressive Retrieval with Guaranteed Error Control*.
- **[SC'25]**: Yanliang Li, Wenbo Li, Qian Gong, Qing Liu, Norbert Podhorszki, Scott Klasky, Xin Liang, and Jieyang Chen. *HP-MDR: High-performance and Portable Data Refactoring and Progressive Retrieval with Advanced GPUs*.
- **[SC'24]**: Xuan Wu, Qian Gong, Jieyang Chen, Qing Liu, Norbert Podhorszki, Xin Liang, and Scott Klasky. *Error-controlled Progressive Retrieval of Scientific Data under Derivable Quantities of Interest*.

ProDM also integrates other existing progressive approaches from researchers and engineers:
- **[TVCG'23]**: Victor AP Magri, Peter Lindstrom. *A general framework for progressive data compression and retrieval*.
- **[SC'21]**: Xin Liang, Qian Gong, Jieyang Chen, Ben Whitney, Lipeng Wan, Qing Liu, David Pugmire, Rick Archibald, Norbert Podhorszki, and Scott Klasky. *Error-controlled, Progressive, and Adaptable Retrieval of Scientific Data with Multilevel Decomposition*.

## Installation

One-command compilation using `build_script.sh`. It will automatically build the ProDM library and its dependencies.

```bash
git clone https://github.com/Linus-Li-1037/ProDM.git
cd ProDM
git switch Predict
sh build_script.sh
```

### Examples

**Example: Hurricane ISABEL dataset**

Hurricane ISABEL data can be downloaded from [SDRBench](https://sdrbench.github.io/). Double-precision versions of `VelocityX.dat`, `VelocityY.dat`, and `VelocityZ.dat` are converted from `Uf48.bin.f32`, `Vf48.bin.f32`, and `Wf48.bin.f32`.

The float versions (`Uf48.bin.f32`, `Vf48.bin.f32`, `Wf48.bin.f32`) are already provided under **`Hurricane/data`** (with the names `VelocityX.dat`, `VelocityY.dat`, and `VelocityZ.dat`) and can be tested directly.

Expected directory layout:

```
Hurricane_d64
├── data
│   ├── VelocityX.dat
│   ├── VelocityY.dat
│   └── VelocityZ.dat
└── refactor
    ├── VelocityX_refactored
    ├── VelocityY_refactored
    └── VelocityZ_refactored
```

**Refactoring and Progressive Retrieval with Multilevel Decomposition [SC'21]**
```bash
cd build
mkdir refactored
# Refactor
# ./test/test_mdr_refactor $data_file $refactored_dict $target_level $num_bitplanes $num_dimensions [dimensions] $dtype
./test/test_mdr_refactor ../Hurricane/data/VelocityX.dat refactored 4 30 3 100 500 500 0 -f
# Retrieval
./test/test_mdr_reconstructor ../Hurricane/data/VelocityX.dat refactored 3 0.01 0.001 0.0001 0 -f
```

**Refactoring and Progressive Retrieval with Iterative Compression [TVCG'23]**
```bash
cd build
mkdir refactored
# Refactor
# ./test/test_pdr_delta_refactor $data_file $refactored_dict $num_dim $dim0 .. $dimn -[dataType: f/d] [Approximator: Dummy-0, MGARD-1, SZ2-2, SZ3-3, HPEZ-4]
./test/test_pdr_delta_refactor ../Hurricane/data/VelocityX.dat refactored 3 100 500 500 -f 3
# Retrieval
# ./test/test_pdr_delta_reconstructor $data_file $refactored_dict num_tolerance tolerance1 ... toleranceN -[dataType: f/d] [Approximator: Dummy-0, MGARD-1, SZ2-2, SZ3-3, HPEZ-4]
./test/test_pdr_delta_reconstructor ../Hurricane/data/VelocityX.dat refactored 3 0.05 0.005 0.0005 -f 3
```

**Progressive Retrieval with QoI error control [SC'24]**
The following steps demonstrate how to test Hurricane ISABEL using `V_total` as the targeted QoI. If the confidential GE data is available, please check the codes in `app/GE` to reproduce the results in the SC'24 paper.  

First convert float data to double for testing and create the refactor directories:
```bash
cd Hurricane
python float2double.py data/VelocityX.dat
python float2double.py data/VelocityY.dat
python float2double.py data/VelocityZ.dat

mkdir Hurricane/refactor
mkdir Hurricane/refactor/VelocityX_refactored
mkdir Hurricane/refactor/VelocityY_refactored
mkdir Hurricane/refactor/VelocityZ_refactored
cd ..
```

Then perform refactoring and retrieval with weighted=0 to disable weighted bitplane encoding (the technique introduced in the HPDC'26 paper):
```bash
cd build
# Refactor
# ./test/refactor_d64 $approximator $weighted=0 $dataset_name $path_to_dataset
./test/refactor_d64 4 0 Hurricane ../Hurricane
# Retrieval
# ./test/qoi_{$target_QoI}_d64 $approximator $weighted $eb $path_to_dataset
./test/qoi_Vtot_d64 4 0 1 $eb ../Hurricane
```

**QoI-based Refactoring and Progressive Retrieval (QProR) [HPDC'26]**

Precision data refactoring using approximators:
```bash
cd build
mkdir refactored
# Refactor
# ./test/test_pdr_refactor $data_file $refactored_dict $num_bitplanes $num_dim $dim0 .. $dimn -[dataType: f/d] [Approximator: Dummy-0, MGARD-1, SZ2-2, SZ3-3, HPEZ-4]
./test/test_pdr_refactor ../Hurricane/data/VelocityX.dat refactored 30 3 100 500 500 -d 4
# Retrieval
# ./test/test_pdr_reconstructor $data_file $refactored_dict num_tolerance tolerance1 ... toleranceN -[dataType: f/d] [Approximator: Dummy-0, MGARD-1, SZ2-2, SZ3-3, HPEZ-4]
./test/test_pdr_reconstructor ../Hurricane/data/VelocityX.dat refactored refactored 3 0.01 0.001 0.0001 -d 4
```
QoI-based refactoring and progressive retrieval with weighted bitplanes (setting weighted=1 to enable weighted bitplane encoding):

```bash
cd build
# Refactor
# ./test/refactor_d64 $approximator $weighted $dataset_id $path_to_dataset $max_weight_v $block_size $approximator_eb
./test/refactor_d64 4 1 Hurricane ../Hurricane 7 4 0.001
# Retrieval
# ./test/qoi_{$target_QoI}_d64 $approximator $weighted $dataset_id $eb $path_to_dataset
./test/qoi_Vtot_d64 4 1 1 $eb ../Hurricane
```
Please refer to  `artifacts/HPDC-26/Appendix.pdf` for artifact description and evaluation.


**Progressive retrieval with Adaptive Interpolation and Coefficient Decomposition (ProAICD) [SC'26]**

```bash
cd build
# Refactor
# ./test/two_modes_refactor data_file output_path -[dataType: f/d] target_level num_bitplanes num_dims dim1 dim2 ... dimn \
#   -[encoder_option: Nega/XOR/PerBit] -[prior_mode: eb(default)/PSNR] -[CP_or_not: CP/no_CP] (coeff_interp_direction, default tune)
./test/two_modes_refactor ../Hurricane/data/VelocityX.dat refactored -d 4 60 100 500 500 -Nega -eb -CP 

# Retrieval
# ./test/two_modes_reconstructor data_file refactored_path -[dataType: f/d] num_of_tolerance tol1 tol2 ... toln \
#   -[encoder_option: Nega/XOR/PerBit] -[interpreter_option: Greedy/DP/BFS] -[CP_or_not: CP/no_CP] [Optional: Reconstructed data path]
./test/two_modes_reconstructor ../Hurricane/data/VelocityX.dat refactored -d 3 0.01 0.001 0.0001 -Nega -DP -CP
```

Please follow  `artifacts/SC-26/evaluation.ipynb` to reproduce the results in the paper.

## Acknoledgment
This project is partially supported by NSF projects under OAC-2628470, OAC-2628472, OAC-2144403, OAC-2311757, OAC-2311758, and DOE RAPIDS-3 SciDAC and Sirius-2 projects. This work used computing resources from Oak Ridge Leadership Computing Facilities (OLCF) and the NSF Advanced Cyberinfrastructure Coordination Ecosystem: Services & Support (ACCESS) program.

## Q&A
Please address your questions to xin.liang@oregonstate.edu with subject title ProDM.
