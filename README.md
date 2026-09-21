# ProDM: A Progressive Data Management Framework for Exascale Science

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

**Prerequisites:** a C++17 compiler, CMake >= 3.18, and libzstd (e.g., `brew install zstd` or `apt install libzstd-dev`). MPI is optional (it is used only by the parallel artifact tools, which are skipped if it is absent), and the example workflow additionally uses python3 with numpy.

One-command compilation using `build_script.sh`. It will automatically build the ProDM library and its dependencies (SZ2, SZ3, QoZ/HPEZ, and MGARD under `external/`), and enables all of them. Compilers default to the system ones and can be overridden, e.g. `CC=gcc-16 CXX=g++-16 sh build_script.sh`.

```bash
git clone https://github.com/lxAltria/ProDM.git
cd ProDM
sh build_script.sh
```

Alternatively, a plain `cmake .. && make` in a build directory produces the dependency-free core (multilevel refactoring, bitplane encoding, and error control: the `mdr` and `proaicd` pipelines; `pdr` and `pdr-delta` need at least one of the SZ3 and HPEZ approximators below). The compressor-based approximators are opt-in CMake options: `-DPRODM_WITH_SZ2=ON`, `-DPRODM_WITH_SZ3=ON`, `-DPRODM_WITH_HPEZ=ON`, and `-DPRODM_WITH_MGARD=ON` (each requires the corresponding library under `external/`, see `build_script.sh`). The command line tools `prodm_refactor` and `prodm_retrieve` (in `test/`) always build, with the approximators that are enabled; the GE application tools under `app/GE` require all four. The evaluation drivers of the published papers live under `artifacts/` (see `artifacts/README.md`) and are built with `-DPRODM_BUILD_ARTIFACTS=ON`, which `build_script.sh` sets.

### Namespaces

All library code lives under the umbrella namespace `ProDM`, organized by what the code is for rather than by which paper introduced it:

- `ProDM` holds the machinery shared by both pipelines: bitplane encoders, level compressors, error control (interfaces, error collectors, the linear estimator `LinearMaxErrorEstimator`, size interpreters), retrievers, writers, and the utilities (`ProDM/Utils`, including the file helpers `readfile`, `writefile`, `print_statistics`).
- `ProDM::MDR` holds the multilevel pipeline (SC'21, SC'26): decomposers, interleavers, tuners, refactors and reconstructors, plus the estimators whose constants come from the multilevel bases (orthogonal basis, cubic interpolation, L2 and s-norm).
- `ProDM::PDR` holds the approximation-based pipeline (TVCG'23, HPDC'26): approximators, refactors, reconstructors.
- `ProDM::MGARDx` holds the in-house multilevel decomposition internals.
- `ProDM::QoI::Hand` holds the QoI registry of the unified command line (`ProDM/QoI/HandQoI.hpp`: QoI values, the hand-derived error bounds of SC'24, the weight recipes of HPDC'26 and the error-bound descents); `ProDM::CLI` the option parsing and file conventions shared by the two tools (`ProDM/Utils/CLI.hpp`).
- `ProDM::Legacy` holds code kept only to reproduce prior papers: the GE synthesizer recipes (`ProDM/App/GE`), `WeightReconstructor` and `QoIRefactor` (`ProDM/Legacy`).

`ProDM/Namespace.hpp` declares the aliases `MDR` and `PDR` at global scope, so `MDR::ComposedRefactor` or `using namespace MDR;` keep compiling; shared components are spelled `ProDM::NegaBinaryBPEncoder`, `ProDM::AdaptiveLevelCompressor` and so on (the former `MDR::` spelling of these no longer compiles, nor does the former `MGARD::` namespace). New headers reopen a namespace with the nested form `namespace ProDM::MDR { ... }` after including `ProDM/Namespace.hpp`.

### Examples

**Example: Hurricane ISABEL dataset**

Hurricane ISABEL data can be downloaded from [SDRBench](https://sdrbench.github.io/). The velocity variables used in the examples are renamed from `Uf48.bin.f32`, `Vf48.bin.f32`, and `Wf48.bin.f32`: single-precision data carries the `.dat.f32` suffix (e.g., `VelocityX.dat.f32`), and the double-precision `.dat` counterparts are derived from them with `float2double.py`. 

```bash
cd example
mkdir -p data
curl -LO https://g-d0cd3f.fd635.8443.data.globus.org/raw-data/Hurricane-ISABEL/SDRBENCH-Hurricane-ISABEL-100x500x500.tar.gz
tar -xvf SDRBENCH-Hurricane-ISABEL-100x500x500.tar.gz
cp 100x500x500/Uf48.bin.f32 data/VelocityX.dat.f32
cp 100x500x500/Vf48.bin.f32 data/VelocityY.dat.f32
cp 100x500x500/Wf48.bin.f32 data/VelocityZ.dat.f32
```

or

```bash
cd example
sh download_data.sh
```

The copied `.dat.f32` files are single-precision (float32) and can be tested directly with the `-f` data type option. 

Expected directory layout (the `.dat` files appear after the float-to-double conversion below):

```
example
├── data
│   ├── VelocityX.dat.f32
│   ├── VelocityX.dat
│   ├── VelocityY.dat.f32
│   ├── VelocityY.dat
│   ├── VelocityZ.dat.f32
│   └── VelocityZ.dat
└── refactor
    ├── VelocityX_refactored
    ├── VelocityY_refactored
    └── VelocityZ_refactored
```

The entire example can be executed by `test_script.sh` after data preparation, and the results are stored in the `*.log` files:

```bash
cd example
sh test_script.sh
```

The following demonstrates a step-by-step breakdown.

**Unified command line: `prodm_refactor` and `prodm_retrieve`**

The two tools `prodm_refactor` and `prodm_retrieve` (sources in `test/`, built under `build/test`) are the common entry
points of all pipelines: `--method` selects the pipeline (`mdr` for multilevel decomposition [SC'21], `proaicd` for
adaptive interpolation with coefficient decomposition [SC'26], `pdr` for approximation-based refactoring [TVCG'23,
HPDC'26], `pdr-delta` for residual snapshots), `--approximator` the compressor behind `pdr` / `pdr-delta`
(`sz3|hpez|ge`, subject to the `PRODM_WITH_*` options), and the remaining options the parameters of
the method. Variables are read from `<data_dir>/<var><suffix>` (`.dat` for `--dtype d`, `.dat.f32` for `--dtype f`)
and each is refactored into `<refactor_dir>/<var>_refactored/`. Without `--qoi`, `prodm_retrieve` retrieves each
variable to tolerances relative to its value range and prints the error statistics; with `--qoi`, the variables of
the QoI are retrieved together under a QoI tolerance.

The QoIs known to `--qoi` and `--weights` are `Vtot`, `Vtot2`, `T`, `C`, `Mach`, `PT` and `mu` over the variables
`VelocityX/Y/Z`, `Pressure`, `Density` (the GE set; `--weights hand:GE` weights both the velocity and the
thermodynamic group, `--max-weight 4,3` gives one maximum weight per group). When the three velocities are refactored
together with `--method pdr`, the nonzero-velocity mask of the SC'24 tools is written to `<refactor_dir>/mask.bin` and
used at retrieval (`--mask none` disables it); the weights and the mask are features of the `pdr` pipeline only. `--joint-range` initializes the per-variable bounds from the joint value range of the QoI's
variables instead of each variable's own range (the convention of the `V_total` tools). The `ge` approximator expects
the GE layout (`<root>/data`, `<root>/refactor`, `<root>/block_sizes.dat`). Run either tool without arguments for the
full option list. The sections below walk through each pipeline with these two tools.


**Refactoring and Progressive Retrieval with Multilevel Decomposition [SC'21]**
```bash
cd build
# Refactor: multilevel decomposition (target level 4) with 30 bitplanes, NegaBinary encoding
./test/prodm_refactor ../example/data refactored --vars VelocityX --dims 100 500 500 --dtype f --method mdr --target-level 4 --bitplanes 30
# Retrieval: tolerances relative to the value range
./test/prodm_retrieve ../example/data refactored --vars VelocityX --dtype f --method mdr --tolerance 0.01 0.001 0.0001
```

**Refactoring and Progressive Retrieval with Iterative Compression [TVCG'23]**
```bash
cd build
# Refactor: residual snapshots on the SZ3 approximator (sz3|hpez|ge)
./test/prodm_refactor ../example/data refactored --vars VelocityX --dims 100 500 500 --dtype f --method pdr-delta --approximator sz3
# Retrieval
./test/prodm_retrieve ../example/data refactored --vars VelocityX --dtype f --method pdr-delta --approximator sz3 --tolerance 0.05 0.005 0.0005
```

**Progressive Retrieval with QoI error control [SC'24]**

The following steps demonstrate how to test Hurricane ISABEL using `V_total` as the targeted QoI. If the confidential GE data is available, please check the codes in `app/GE` (and `artifacts/SC-24`, `artifacts/HPDC-26`) to reproduce the results of the SC'24 and HPDC'26 papers.  

First convert float data to double for testing:
```bash
cd example
python float2double.py data/VelocityX.dat.f32
python float2double.py data/VelocityY.dat.f32
python float2double.py data/VelocityZ.dat.f32
```

Then refactor the three velocities without weights and retrieve them under a QoI tolerance (relative to the value range of `V_total`); the hand-derived estimator bounds the QoI error from the per-variable bounds, and the coordinate descent tightens the bounds until the estimate meets the tolerance:
```bash
cd build
# Refactor (the nonzero-velocity mask is written to ../example/refactor/mask.bin)
./test/prodm_refactor ../example/data ../example/refactor --vars VelocityX,VelocityY,VelocityZ --dims 100 500 500 --dtype d --method pdr --approximator hpez
# Retrieval (--vars defaults to the variables of the QoI)
./test/prodm_retrieve ../example/data ../example/refactor --dtype d --method pdr --approximator hpez --qoi Vtot --estimator hand --inverse coordinate --joint-range --tolerance 0.01
```

**QoI-based Refactoring and Progressive Retrieval (QProR) [HPDC'26]**

Precision data refactoring using approximators:
```bash
cd build
# Refactor: approximation-based refactoring on HPEZ (approximator bound 0.001, 30 bitplanes)
./test/prodm_refactor ../example/data refactored --vars VelocityX --dims 100 500 500 --dtype d --method pdr --approximator hpez --bitplanes 30
# Retrieval
./test/prodm_retrieve ../example/data refactored --vars VelocityX --dtype d --method pdr --approximator hpez --tolerance 0.01 0.001 0.0001
```
QoI-based refactoring and progressive retrieval with weighted bitplanes (`--weights hand:<qoi>` derives per-point weights from the QoI; `--max-weight` and `--block-size` are the weighting parameters, `--eb` the approximator bound):

```bash
cd build
# Refactor
./test/prodm_refactor ../example/data ../example/refactor --vars VelocityX,VelocityY,VelocityZ --dims 100 500 500 --dtype d --method pdr --approximator hpez --weights hand:Vtot --eb 0.001 --max-weight 7 --block-size 4
# Retrieval (the stored weights are detected; the coordinate descent becomes the proportional update of QProR)
./test/prodm_retrieve ../example/data ../example/refactor --dtype d --method pdr --approximator hpez --qoi Vtot --estimator hand --inverse coordinate --joint-range --tolerance 0.01
```
Please refer to `artifacts/HPDC-26/README.md` for the artifact description and evaluation instructions.


**Progressive retrieval with Adaptive Interpolation and Coefficient Decomposition (ProAICD) [SC'26]**

```bash
cd build
# Refactor: encoder nega|xor|perbit, prior eb|psnr, --cp enables coefficient decomposition
./test/prodm_refactor ../example/data refactored --vars VelocityX --dims 100 500 500 --dtype d --method proaicd --target-level 4 --bitplanes 60 --encoder nega --prior eb --cp
# Retrieval: interpreter greedy|dp; pass the same --encoder and --cp as the refactor
./test/prodm_retrieve ../example/data refactored --vars VelocityX --dtype d --method proaicd --encoder nega --interpreter dp --cp --tolerance 0.01 0.001 0.0001
```

Please refer to `artifacts/SC-26/README.md` and `artifacts/SC-26/ablation_steps.sh` to reproduce the results in the paper.

## Acknowledgment
This project is partially supported by NSF projects under OAC-2628470, OAC-2628472, OAC-2144403, OAC-2311757, OAC-2311758, and DOE RAPIDS-3 SciDAC and Sirius-2 projects. This work used computing resources from Oak Ridge Leadership Computing Facilities (OLCF) and the NSF Advanced Cyberinfrastructure Coordination Ecosystem: Services & Support (ACCESS) program. This work used Claude Code for code refactoring and review. 

## Q&A
Please address your questions to xin.liang@oregonstate.edu with subject title ProDM.
