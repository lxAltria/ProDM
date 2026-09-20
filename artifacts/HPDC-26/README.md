# HPDC'26: QProR, an Efficient Framework for Quantity-of-Interest Based Progressive Retrieval with Guaranteed Error Control

Wenbo Li, Qian Gong, Xuan Wu, Jieyang Chen, Qing Liu, Xubin He, Norbert Podhorszki, Scott Klasky, Xin Liang. HPDC'26.

References: the paper's full AD/AE are located in https://github.com/Linus-Li-1037/ProDM/tree/QPro (Appendix.pdf) and https://github.com/Linus-Li-1037/qoi-control/tree/qpro (code repository). This folder documents the components integrated into ProDM.

## What this folder holds

- `parallel/para_refactor`, `parallel/para_Vtot`: the parallel drivers (MPI; float data; approximator selectable; weighted or unweighted encoding; uniform or coordinate descent), built with `-DPRODM_BUILD_ARTIFACTS=ON`.
- `parallel/paraGE_{Vtot,T,C,Mach,PT,mu}`, `parallel/para_refactor_d64`: the older double-precision GE drivers, kept as sources only; superseded by the two above.
- `structured_steps.sh`: the QProR-N / QProR-NW / QProR ladders on a structured dataset (see below).

The sequential tools themselves are the GE application in `app/GE`: `refactor_d64` and `qoi_{Vtot,T,C,Mach,PT,mu}_d64`.

## Method names and tool arguments

| method | refactoring | retrieval |
|---|---|---|
| QProR-N: approximation + bitplane encoding, uniform tightening | `refactor_d64 A 0 <Dataset> <Path>` | `qoi_<QoI>_d64 A 0 0 <eb> <Path>` |
| QProR-NW: + weighted bitplane encoding | `refactor_d64 A 1 <Dataset> <Path> <w> ...` | `qoi_<QoI>_d64 A 1 0 <eb> <Path>` |
| QProR: + coordinate error-bound descent | same representation as QProR-NW | `qoi_<QoI>_d64 A 1 1 <eb> <Path>` |

`A` is the approximator: 3 = the GE approximator (GE data), 1 = SZ3-based approximation (structured data); 4 = HPEZ, 2 = PMGARD, 5 = SZ2, 6 = MGARD, 0 = Dummy are also available. The weighted refactoring takes `<max_weight_V> <max_weight_T> <approximator_eb>` for GE and `<max_weight> <block_size> <approximator_eb>` for structured data; the paper used 4 3 0.001 on GE and 7 4 0.001 elsewhere. The third retrieval argument selects the descent (0 uniform, 1 coordinate); no source editing is needed, unlike the appendix's instructions, which predate this option.

## Reproducible with this repository

**GE, all six QoIs, QProR-N / QProR-NW / QProR** (Figures 4 to 6, the QProR rows): build with `sh build_script.sh`, arrange the reordered GE data as `<GE>/data/{VelocityX,VelocityY,VelocityZ,Pressure,Density}.dat` with `block_sizes.dat` in `<GE>/`, and run

    DATA_DIR=<GE> sh app/GE/GE_steps.sh

which writes `build/Result/GE_result_3D.txt` (Methods 4, 5, 6 are QProR-N, QProR-NW, QProR; the twelve error bounds 0.1 to 5e-7 relative to the QoI range). Individual points: e.g. `./build/app/GE/qoi_Vtot_d64 3 1 1 1e-3 <GE>`.

**Structured datasets, V_total** (Figure 7, the QProR rows): arrange a dataset as `<D>/data/{VelocityX,VelocityY,VelocityZ}.dat` (double precision; Hurricane 100x500x500, NYX 512^3, SCALE 98x1200x1200, S3D 500^3 are known to the tools by name) and run

    DATASET=Hurricane DATA_DIR=<D> sh artifacts/HPDC-26/structured_steps.sh

for the three methods over the same twelve error bounds (`APPROX`, `MAX_WEIGHT`, `BLOCK_SIZE` override the paper's 1, 7, 4). Results go to `build/Result/<Dataset>_result.txt`.

**Parallel runs**: `mpirun -n <cores> ./build/artifacts/HPDC-26/para_refactor ...` followed by `para_Vtot [Approximator] [weighted] [decrease_method] [eb] [input path] [output path]`.

## Requires the qoi-control repository

The baselines PMGARD-HB and PSZ3-delta (Figures 6 to 8) are the `refactor_data`, `halfing_Vtot` and `halfing_Vtot_sz3delta` tools of https://github.com/Linus-Li-1037/qoi-control/tree/qpro. `prodm_retrieve --qoi Vtot` does run the same QoI loop over ProDM's multilevel (`--method mdr`) and delta (`--method pdr-delta`) pipelines, but these are not the qoi-control tools the figures show (different stream layout, no velocity mask in the multilevel path) and are not a like-for-like substitute. IPComp is a separate library.

## Datasets

GE-small and GE-large (unstructured) are private data that could not be shared without permission; Hurricane, NYX, SCALE and S3D are public (for S3D the last three fields of the last snapshot are the velocities). Hurricane, NYX and SCALE were converted from single to double precision (`example/float2double.py`).

## Differences from the published numbers

The hand-derived bounds (`include/ProDM/Legacy/QoIUtils.hpp`) return +infinity where a bound's precondition fails, where the submission-time code returned 0 and silently accepted the point. Such a point is now the worst point of the pass, so the retrieval keeps tightening the data bounds until its bound becomes finite and meets the tolerance (the uniform descent divides every bound by 1.5, and the coordinate descent also tightens every variable by the full factor while the estimate is infinite, since an infinite estimate gives it no direction). On the GE ladders this changes 7 of the 216 retrievals, at the error bounds 0.1 and 0.05 of Mach, PT and mu under QProR-NW and QProR, by at most 0.56 bits, with no tolerance violation before or after.
