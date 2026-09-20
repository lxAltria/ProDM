# SC'24: Error-controlled Progressive Retrieval of Scientific Data under Derivable Quantities of Interest

Xuan Wu, Qian Gong, Jieyang Chen, Qing Liu, Norbert Podhorszki, Xin Liang, Scott Klasky. SC'24.

References: the paper's implementation and artifact description are in the qoi-control repository
(https://github.com/xuanwu02/qoi-control). This folder documents the components integrated into ProDM.

## What this folder holds

`parallel/`: the SC'24 parallel GE drivers, MPI with ADIOS2 I/O (`.bp` streams). Each pair is one method of the paper
on V_total:

| method | refactoring | retrieval |
|---|---|---|
| PMGARD-HB: multilevel decomposition (hierarchical basis) with bitplane encoding | `refactorGE` | `testGE_VTOT` |
| SZ3 multi-precision: independent SZ3 snapshots at decreasing error bounds | `refactorGE_SZ3` | `testGE_VTOT_SZ3` |
| PSZ3-delta: SZ3 snapshots of the residuals | `refactorGE_SZ3_delta` | `testGE_VTOT_SZ3_delta` |

They build with `-DPRODM_BUILD_ARTIFACTS=ON` when MPI, ADIOS2 and SZ3 are available; build_script.sh does not build
ADIOS2 (see its commented block), so they are skipped by default.

## Reproducible with this repository

**PMGARD-HB on GE, sequential.** The GE tools in `app/GE` (the GE application) run the multilevel pipeline inside the
SC'24 uniform-tightening loop with the hand-derived bounds when the approximator id is 2. Arrange the reordered GE
data as `<GE>/data/{VelocityX,VelocityY,VelocityZ,Pressure,Density}.dat` and run

    ./build/app/GE/refactor_d64 2 0 GE <GE>
    for eb in 0.1 0.05 0.01 0.005 0.001 0.0005 0.0001 0.00005 0.00001 0.000005 0.000001 0.0000005; do
        ./build/app/GE/qoi_Vtot_d64 2 0 0 $eb <GE>
    done

Each retrieval prints the requested and estimated QoI errors, the actual error and the bitrate. The same
representation serves the other five GE QoIs through `qoi_{T,C,Mach,PT,mu}_d64 2 0 0 <eb> <GE>`. The refactoring
recipe is one-dimensional (`refactor_velocities_1D_PMGARD_BP`), as in the paper; `refactor_d64 2 0 Hurricane|NYX`
gives the 3D variant on those datasets.

`app/GE/GE_steps.sh` runs the HPDC'26 methods (QProR-N, QProR-NW, QProR) on the same data; see `artifacts/HPDC-26`.

## Requires the qoi-control repository

The sequential SZ3 multi-precision and PSZ3-delta retrievals are in qoi-control only; ProDM has their refactoring
recipes (`refactor_GE_SZ3`, `refactor_GE_SZ3_delta` in the GE synthesizer) and the parallel ADIOS2 drivers above, but
no sequential QoI-controlled retrieval over SZ3 snapshots or residuals.

## Differences from the published numbers

The hand-derived bounds (`include/ProDM/Legacy/QoIUtils.hpp`) return +infinity where a bound's precondition fails
(a divisor whose interval contains zero), where the original code returned 0. This changes the forward estimate: a
point at which a precondition fails used to be scored 0 and ignored, and is now the worst point of the pass, so the
loop tightens the data bounds until its bound becomes finite and meets the tolerance, under any descent. The SC'24
results are unaffected because the V_total chain (squares and a square root) has no failing branch. The QoIs with
divisions (T, C, Mach, PT, mu) can change: in the GE reproduction, the unweighted uniform runs (QProR-N) kept all
their rows because the preconditions held at the initial bounds, while the weighted runs changed on Mach, PT and mu
at the two loosest tolerances (see `artifacts/HPDC-26`).
