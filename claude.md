# Project: ProDM -- A Progressive Data Refactoring and Retrieval Framework for Scientific Applications

Develop a progressive framework for scientific data management. Include both sequential code for serial execution and parallel code for distributed execution with MPI.

## Build

- `mkdir build && cd build && cmake .. && make -j` — core build (multilevel decomposition, bitplane encoding, error control; requires only a C++17 compiler and system zstd)
- `sh build_script.sh` — full build: fetches and builds the optional compressor dependencies (SZ2, SZ3, QoZ/HPEZ, MGARD) under `external/` and configures with all `PRODM_WITH_*` options ON
- Optional CMake options (default OFF): `PRODM_WITH_SZ2`, `PRODM_WITH_SZ3`, `PRODM_WITH_HPEZ`, `PRODM_WITH_MGARD` (approximators); `PRODM_WITH_MPI` (default ON, skipped gracefully if MPI is absent)
- Compilers default to the system `cc`/`c++`; override with `CC=gcc-16 CXX=g++-16 sh build_script.sh`. MGARD exposes a C++ API, so `libmgard` must be built with the same C++ toolchain as ProDM.

### Dependencies

- CMake >= 3.13
- C++17 compiler
- zstd (system package, e.g. `brew install zstd` / `apt install libzstd-dev`)
- MPI (parallel tools only; auto-detected, skipped if absent)
- python3 with numpy (for `example/float2double.py` used by the example workflow)

## Structure

- `app`           — code for specific applications (GE); not built by default
- `artifacts`     — paper artifacts (Appendix, evaluation notebooks, reproducibility scripts)
- `example`       — example dataset workflow (Hurricane ISABEL): `download_data.sh`, `test_script.sh`, `float2double.py`
- `include/ProDM` — header-only library, organized by pipeline stage (Decomposer/{MultiLevel,Approximation}, Encoder, Compressor, Writer, Retriever, ErrorControl/{Collector,Estimator,SizeInterpreter}, Refactor/{MDR,PDR}, Reconstructor/{MDR,PDR}, Utils, App/GE)
- `prl_src`       — parallel (MPI) drivers
- `test`          — sequential drivers for evaluation

## Code Conventions

### Memory and data layout
- All arrays use row-major layout.
- Single-precision example data uses the `.dat.f32` suffix; double-precision uses `.dat` (derived via `float2double.py`).

### General C++ style
- C++17 only. Do not use C++20/23 features.
- Use `const` and pass by reference wherever possible.

## Debugging and Testing

### Evaluation
- Run `sh test_script.sh` from `example/` (after data download) to exercise all workflows; results land in `example/*.log`.
- Check the python notebook in `artifacts/SC-26/evaluation.ipynb` for paper results.

## Things NOT to Do
- Never modify `external` and the data under `example/data`.
- Do not modify `app` for now.
