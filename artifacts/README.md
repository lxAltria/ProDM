# Artifacts

Evaluation drivers and instructions for the published papers. The library and the application tools (`include/`,
`test/`, `app/`) are what ProDM ships; the folders here are kept so that each paper's experiments can be repeated.
Build them with `-DPRODM_BUILD_ARTIFACTS=ON` (build_script.sh does), under the same approximator options as the
library; MPI tools need `PRODM_WITH_MPI=ON` (default) and an MPI installation.

| folder | paper | contents |
|---|---|---|
| `SC-24` | Wu et al., SC'24: error-controlled progressive retrieval under derivable QoIs (GE, V_total) | parallel GE drivers with ADIOS2 I/O; pointers to the GE tools in `app/GE` |
| `HPDC-26` | Li et al., HPDC'26: QProR | `structured_steps.sh` (the QProR ladders on a structured dataset) and the parallel weighted refactoring and QoI retrieval drivers (current float tools and the older double-precision GE tools); the artifact-evaluation instructions of the appendix |
| `SC-26` | Li et al., SC'26: adaptive interpolation and coefficient decomposition (proaicd) | `ablation_steps.sh` (dataset preparation, the error-bound and PSNR ladders, Tables 3 to 4, the Figure 12 reconstructions) and the parallel drivers of the multilevel, proaicd and PSZ3-delta pipelines |

SC'25 (HP-MDR) is evaluated in the official MGARD repository and has no folder here.
