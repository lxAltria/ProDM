# ProDM: Progressive Data Refactoring and Retrieval Framework

This repository hosts the code for two related papers:

1. **[HPDC'26]**: *QProR: An Efficient Framework for Quantity-of-Interest Based Progressive Retrieval with Guaranteed Error Control*
2. **[SC'26]** *Improving Progressive Compression with Adaptive Interpolation and Coefficient Decomposition*

---

## Authors

**Major authors:**
Dr. Xin Liang (Oregon State University), Dr. Qing Liu (NJIT), Dr. Xubin He (Temple)

**Other contributors:**
Xuan Wu (Oregon State University), Qirui Tian (NJIT), Wenbo Li (Oregon State University)

**Collaborators:**
Dr. Scott Klasky (ORNL), Dr. Qian Gong (ORNL), Dr. Jill Zhang (LLNL), Dr. Seung-Hoe Ku (PPPL), Dr. Xiaohua Zhang (LLNL), Dr. Jieyang Chen (UAB), etc.

---

## Paper 1: [HPDC'26] QProR: An Efficient Framework for Quantity-of-Interest Based Progressive Retrieval with Guaranteed Error Control

### Installation

One-command compilation using `build_script.sh`. It will automatically build the ProDM library and its dependencies.

```bash
git clone https://github.com/Linus-Li-1037/ProDM.git
cd ProDM
git switch Predict
sh build_script.sh
```

### Examples

**HPDC'26 Template**

Precision data refactoring using approximators:

```bash
cd build
# Refactor
./test/refactor_d64 $approximator $wbp $dataset $path_to_dataset $max_weight_v $block_size
# Retrieval
./test/qoi_{$target_QoI}_d64 $approximator $wbp $eb $path_to_dataset
```

**Example: Hurricane ISABEL dataset**

Hurricane ISABEL data can be downloaded from [SDRBench](https://sdrbench.github.io/). Double-precision versions of `VelocityX.dat`, `VelocityY.dat`, and `VelocityZ.dat` are extended from `Uf48.bin.f32`, `Vf48.bin.f32`, and `Wf48.bin.f32`.

The float versions (`Uf48.bin.f32`, `Vf48.bin.f32`, `Wf48.bin.f32`) are already provided under **`Hurricane/data`** and can be tested directly.

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

Commands to test Hurricane ISABEL using `V_total` as the targeted QoI:

```bash
python float2double.py Hurricane/data/VelocityX.dat
python float2double.py Hurricane/data/VelocityY.dat
python float2double.py Hurricane/data/VelocityZ.dat

mkdir Hurricane/refactor
mkdir Hurricane/refactor/VelocityX_refactored
mkdir Hurricane/refactor/VelocityY_refactored
mkdir Hurricane/refactor/VelocityZ_refactored

cd build
# Refactor
./test/refactor_d64 4 1 Hurricane ../Hurricane 7 4
# Retrieval
./test/qoi_Vtot_d64 4 1 1 $eb ../Hurricane
```

### Artifacts

Please check `Appendix.pdf` for a detailed description.

---

## Paper 2: [SC'26] Improving Progressive Compression with Adaptive Interpolation and Coefficient Decomposition

### Installation

One-command compilation using `build_script.sh`. It will automatically build the ProDM library and its dependencies.

```bash
git clone https://github.com/Linus-Li-1037/ProDM.git
cd ProDM
git switch Predict
sh build_script.sh
```

### Examples

**SC'26 Template**

```bash
cd build
# Refactor
./test/two_modes_refactor data_file -[dataType: f/d] target_level num_bitplanes num_dims dim1 dim2 ... dimn output_path \
  -[encoder_option: Nega/XOR/PerBit] -[prior_mode: eb(default)/PSNR] -[CP_or_not: CP/no_CP] (coeff_interp_direction, default tune)

# Retrieval
./test/two_modes_reconstructor data_file -[dataType: f/d] num_of_tolerance tol1 tol2 ... toln refactored_path \
  -[encoder_option: Nega/XOR/PerBit] -[interpreter_option: Greedy/DP/BFS] -[CP_or_not: CP/no_CP] [Optional: Reconstructed data path]
```

### Artifacts

Please use `evaluation.py` to reproduce the results in the paper.