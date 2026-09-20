# SC'26: Improving Progressive Compression with Adaptive Interpolation and Coefficient Decomposition

Wenbo Li, Xuan Wu, Qian Gong, Pu Jiao, Jieyang Chen, Qing Liu, Norbert Podhorszki, Scott Klasky, Xin Liang. SC'26.

References: the artifact of the paper (job scripts, evaluation notebook, appendix) is located at https://github.com/Linus-Li-1037/ProDM/tree/Predict (Zenodo record 22102905). This folder
documents the components integrated into ProDM and how the paper's sequential results are reproduced with them.

## What this folder holds

- `ablation_steps.sh`: the reproduction script. It prepares a dataset, runs the error-bound and PSNR ladders of the
  proposed method and of the PMGARD-HB and SZ3-R baselines with `prodm_refactor` / `prodm_retrieve`, summarizes
  them as in Tables 3 and 4, and produces the CESM temperature reconstructions of Figure 12.
- `parallel/`: the MPI drivers of the three pipelines used in the weak-scaling study (Figure 11), built with
  `-DPRODM_BUILD_ARTIFACTS=ON` when MPI is found: `para_proaicd_refactor` / `para_proaicd_reconstructor`,
  `para_mdr_refactor` / `para_mdr_reconstructor`, `para_PSZ3-delta_refactor` / `para_PSZ3-delta_reconstructor`
  (under `build/artifacts/SC-26/`).

## What can be reproduced

| paper element | status |
|---|---|
| Figures 7 to 10 (error-bound and PSNR ladders on S3D, CESM, Miranda, SCALE) | reproduced by `ablation_steps.sh ablation` for the proposed method (AdatInterp, CoeffDecom, both modes), PMGARD-HB and SZ3-R; the IPComp curves are not |
| Tables 3 and 4 (average refactoring time; retrieval time and bitrate at tolerance 1e-4) | `ablation_steps.sh summary`, same three methods |
| Figure 12 (CESM temperature field) | `ablation_steps.sh visualize` writes the PMGARD-HB, SZ3-R and proposed reconstructions; the IPComp panel is not |
| Figure 6 (decomposition and encoding throughput) | not reproducible here: the measurement used `test_decompose` of the external MGARDx build and a `test_mdr_fuse_refactor` timing tool that is not part of ProDM |
| Figure 11 (weak scaling on JHTDB, 1024 cores) | the three ProDM MPI drivers are in `parallel/`; the study additionally needs the external IPComp MPI tools, the 512 GB JHTDB cubes, a 1024-core homogeneous partition and a Globus transfer, and is not scripted here |

IPComp (Magri and Lindstrom, TVCG'23) is an external library that is not included in ProDM; every IPComp row,
curve and panel of the paper is outside this reproduction.

## Methods and tools

Every method is one `prodm_refactor` / `prodm_retrieve` configuration (`D` and `R` are the dataset's `data` and
`refactor` directories, `V` the variable; `D/V.dat` is refactored into `R/V_refactored/`):

| # | method | refactoring | retrieval |
|---|---|---|---|
| 2 | PMGARD-HB, eb mode | `prodm_refactor D R --vars V --dims n1 n2 n3 --dtype d --method mdr --target-level L --bitplanes 60 --encoder perbit` | `prodm_retrieve D R --vars V --dtype d --method mdr --encoder perbit --tolerance t` |
| 3 | PMGARD-HB, PSNR mode | same with `--encoder nega` | same with `--encoder nega` |
| 4 | SZ3-R (residual snapshots on SZ3) | `prodm_refactor D R --vars V --dims n1 n2 n3 --dtype d --method pdr-delta --approximator sz3` | `prodm_retrieve D R --vars V --dtype d --method pdr-delta --approximator sz3 --tolerance t` |
| 5 | AdatInterp, eb mode | `prodm_refactor D R --vars V --dims n1 n2 n3 --dtype d --method proaicd --target-level L --bitplanes 60 --encoder perbit --prior eb` | `prodm_retrieve D R --vars V --dtype d --method proaicd --encoder perbit --interpreter dp --tolerance t` |
| 6 | CoeffDecom, eb mode (proposed) | same with `--cp` | same with `--cp` |
| 7 | AdatInterp, PSNR mode | `--encoder xor --prior psnr` | `--encoder xor --interpreter dp` |
| 8 | CoeffDecom, PSNR mode (proposed) | `--encoder xor --prior psnr --cp` | `--encoder xor --interpreter dp --cp` |

The numbers are those of the paper's scripts (#1 was IPComp). Tolerances are relative to the value range of the
variable; the target level `L` is 3 for CESM and 4 for the other datasets; every bound is retrieved from scratch in
its own process, as in the paper's measurement.

## Data

The four datasets come from [SDRBench](https://sdrbench.github.io/) and are stored as one double-precision file per
variable, `<DATA_ROOT>/<dataset>/data/<var>.dat`, with the refactored streams under `<DATA_ROOT>/<dataset>/refactor/`.
`ablation_steps.sh prepare` downloads and converts a dataset (python3 with numpy; the archives are removed after
conversion):

| dataset | dimensions | variables (fields evaluated in the paper) | raw size |
|---|---|---|---|
| CESM-ATM | 26 x 1800 x 3600 | 33 fields (`CLDICE` ... `Z3`; float32 source) | 41 GB |
| Miranda | 256 x 384 x 384 | density, diffusivity, pressure, velocityx, velocityy, velocityz, viscocity | 2 GB |
| SCALE-LETKF | 98 x 1200 x 1200 | PRES, QC, QG, QI, QR, QS, QV, RH, U, V, W (float32 source) | 12 GB |
| S3D | 500 x 500 x 500 | CH4, CO, CO2, H2O, O2, Temperature, VelocityX, VelocityY, VelocityZ (sliced from the 11-variable field file) | 8 GB |

With the refactored representations of all methods the four datasets take roughly 200 GB.

## Running

Build ProDM with `sh build_script.sh` (SZ3 is required for SZ3-R; `PRODM_WITH_SZ3=ON`). Then, per dataset:

```bash
cd artifacts/SC-26
DATASET=S3D DATA_ROOT=/path/to/datasets bash ablation_steps.sh prepare      # once
DATASET=S3D DATA_ROOT=/path/to/datasets bash ablation_steps.sh ablation     # all variables, 7 methods, 17 bounds
DATASET=S3D DATA_ROOT=/path/to/datasets bash ablation_steps.sh summary      # Tables 3 and 4 for this dataset
DATASET=CESM DATA_ROOT=/path/to/datasets bash ablation_steps.sh visualize   # Figure 12 reconstructions
```

`ablation` writes one result file per variable, `build/Result/SC-26/<dataset>_<var>_ablation.txt` (override with
`RESULT_DIR`), with one line per method and bound:

```
Method#6, CoeffDecom-EB, Refactor time: 12.3s
Method#6, CoeffDecom-EB, ErrorBound=0.0001, Bitrate = 2.87, MSE = ..., PSNR = ..., NRMSE = ..., Max error = ..., Reconstruct time: 1.1s
```

The bitrate-versus-bound and bitrate-versus-PSNR curves of Figures 7 to 10 are these lines per variable; the paper
plots the bounds with mantissa 1 and 5 (all 17 of the ladder). Variables are independent: `VARS="CH4 CO"` restricts a
run, so the fields of a dataset can be spread over processes or scheduler jobs (the paper's runs used one job per
field). A full dataset takes from about 10 minutes (Miranda) to about 40 minutes (CESM) per field-parallel run; the
S3D ladder is the most expensive per field. `EBS="0.1 0.01"` shortens the ladder for a test.

`visualize` (CESM only) reconstructs the `T` field with PMGARD-HB in PSNR mode at 4e-3, SZ3-R at 1e-2 and 1e-3, and the
proposed method in PSNR mode at 2e-3, and writes `PMGARD_T.dat`, `SZ3R_T_low.dat`, `SZ3R_T_high.dat` and `TWO_T.dat`
(26 x 1800 x 3600 doubles) into `<DATA_ROOT>/CESM/refactor/T_refactored/`, with their bitrates and PSNR in
`build/Result/SC-26/CESM_T_visualization.txt`; Figure 12 shows slice 11, rows 1000 to 1200, columns 1480 to 1680.

## MPI drivers (Figure 11)

The weak-scaling study refactored 1024 cubes of the JHTDB `rotstrat4096` temperature (256 x 512 x 512 doubles,
`Temperature_<rank>.d64`) with one rank per cube and retrieved them at the relative tolerance 1e-4 on 128 to 1024
cores; the drivers take the cube directory and the rank as the file index:

```bash
mpirun -n <ranks> build/artifacts/SC-26/para_mdr_refactor <data>/Temperature_ 4 60 3 256 512 512 <refactor>/
mpirun -n <ranks> build/artifacts/SC-26/para_mdr_reconstructor <data>/Temperature_ 1 1e-4 <refactor>/ <output>/
mpirun -n <ranks> build/artifacts/SC-26/para_PSZ3-delta_refactor <data>/Temperature_ <refactor>/ 18 3 256 512 512 -d
mpirun -n <ranks> build/artifacts/SC-26/para_PSZ3-delta_reconstructor <data>/Temperature_ <refactor>/ 1 1e-4 -d <output>/
mpirun -n <ranks> build/artifacts/SC-26/para_proaicd_refactor <data>/Temperature_ -d 4 60 3 256 512 512 <refactor>/ -PerBit -eb -CP
mpirun -n <ranks> build/artifacts/SC-26/para_proaicd_reconstructor <data>/Temperature_ -d 1 1e-4 <refactor>/ -PerBit -DP -CP <output>/
```

Each driver prints `max_elapsed_time` (the slowest rank); the paper adds the Globus transfer time of the retrieved
data. The drivers keep the argument lists of the paper's `para_*` tools (see each usage line) and are not routed
through `prodm_refactor` / `prodm_retrieve`.
