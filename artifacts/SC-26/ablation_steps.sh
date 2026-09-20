#!/bin/bash
# ============================================================================
# SC'26 reproduction with prodm_refactor / prodm_retrieve:
#   error-bound and PSNR ladders of the proposed method (adaptive interpolation,
#   coefficient decomposition) and the PMGARD-HB and SZ3-R baselines on one
#   dataset (Figures 7 to 10, Tables 3 to 4), and the CESM temperature
#   reconstructions of Figure 12. IPComp is an external library and is not run.
#
# Usage:  DATASET=<CESM|Miranda|SCALE|S3D> DATA_ROOT=<dir> bash ablation_steps.sh [stage ...]
# Stages: prepare    download the dataset from SDRBench into ${DATA_ROOT}/${DATASET}/data (python3 + numpy)
#         ablation   run the seven-method ladder on every variable (default stage)
#         summary    average refactoring time, and retrieval time and bitrate at the tolerance 1e-4 (Tables 3, 4)
#         visualize  CESM only: reconstruct the temperature field at the bounds of Figure 12
#
# Layout: ${DATA_ROOT}/${DATASET}/data/<var>.dat (double precision), refactored streams in
#         ${DATA_ROOT}/${DATASET}/refactor/<var>_refactored/, results in ${RESULT_DIR}/${DATASET}_<var>_ablation.txt
# Overrides: VARS="A B" (subset of variables; one variable per process is the natural unit for GNU parallel or
#            a scheduler), EBS="0.1 0.01" (error-bound ladder), DIMS / LEVEL (with DATASET=custom), PRODM_DIR,
#            RESULT_DIR. Requires ProDM built with PRODM_WITH_SZ3=ON (build_script.sh enables it).
# ============================================================================
set -u
script_dir=$(cd "$(dirname "$0")" && pwd)
PRODM_DIR=${PRODM_DIR:-$(cd "${script_dir}/../.." && pwd)}
R=${PRODM_DIR}/build/test/prodm_refactor
Q=${PRODM_DIR}/build/test/prodm_retrieve
DATASET=${DATASET:?set DATASET to CESM, Miranda, SCALE, S3D or custom}
DATA_ROOT=${DATA_ROOT:?set DATA_ROOT to the directory holding <DATASET>/data}
RESULT_DIR=${RESULT_DIR:-${PRODM_DIR}/build/Result/SC-26}
D=${DATA_ROOT}/${DATASET}/data
RF=${DATA_ROOT}/${DATASET}/refactor

case ${DATASET} in
    CESM)    DIMS=${DIMS:-"26 1800 3600"}; LEVEL=${LEVEL:-3}
             VARS=${VARS:-"CLDICE CLDLIQ CLOUD CMFDQR CMFDQ CMFDT CONCLD DCQ DTCOND DTV FICE GCLDLWP ICIMR ICLDIWP ICLDTWP ICWMR OMEGAT OMEGA QC QRL QRS Q RELHUM T UU U VD01 VQ VT VU VV V Z3"} ;;
    Miranda) DIMS=${DIMS:-"256 384 384"}; LEVEL=${LEVEL:-4}
             VARS=${VARS:-"density diffusivity pressure velocityx velocityy velocityz viscocity"} ;;
    SCALE)   DIMS=${DIMS:-"98 1200 1200"}; LEVEL=${LEVEL:-4}
             VARS=${VARS:-"PRES QC QG QI QR QS QV RH U V W"} ;;
    S3D)     DIMS=${DIMS:-"500 500 500"}; LEVEL=${LEVEL:-4}
             VARS=${VARS:-"CH4 CO CO2 H2O O2 Temperature VelocityX VelocityY VelocityZ"} ;;
    custom)  DIMS=${DIMS:?set DIMS for a custom dataset}; LEVEL=${LEVEL:?set LEVEL for a custom dataset}; VARS=${VARS:?set VARS for a custom dataset} ;;
    *) echo "unknown DATASET ${DATASET}"; exit 1 ;;
esac
# the 17 relative error bounds of the paper's ladders
EBS=${EBS:-"0.1 0.05 0.01 0.005 0.001 0.0005 0.0001 0.00005 0.00001 0.000005 0.000001 0.0000005 0.0000001 0.00000005 0.00000001 0.000000005 0.000000001"}

# ---------------------------------------------------------------------------
# method table (the numbers follow the paper's scripts; #1 was IPComp)
# ---------------------------------------------------------------------------
method_name()   { case $1 in 2) echo PMGARD-EB;; 3) echo PMGARD-PSNR;; 4) echo SZ3-R;; 5) echo AdatInterp-EB;; 6) echo CoeffDecom-EB;; 7) echo AdatInterp-PSNR;; 8) echo CoeffDecom-PSNR;; esac; }
refactor_opts() { case $1 in
    2) echo "--method mdr --target-level ${LEVEL} --bitplanes 60 --encoder perbit";;
    3) echo "--method mdr --target-level ${LEVEL} --bitplanes 60 --encoder nega";;
    4) echo "--method pdr-delta --approximator sz3";;
    5) echo "--method proaicd --target-level ${LEVEL} --bitplanes 60 --encoder perbit --prior eb";;
    6) echo "--method proaicd --target-level ${LEVEL} --bitplanes 60 --encoder perbit --prior eb --cp";;
    7) echo "--method proaicd --target-level ${LEVEL} --bitplanes 60 --encoder xor --prior psnr";;
    8) echo "--method proaicd --target-level ${LEVEL} --bitplanes 60 --encoder xor --prior psnr --cp";; esac; }
retrieve_opts() { case $1 in
    2) echo "--method mdr --encoder perbit";;
    3) echo "--method mdr --encoder nega";;
    4) echo "--method pdr-delta --approximator sz3";;
    5) echo "--method proaicd --encoder perbit --interpreter dp";;
    6) echo "--method proaicd --encoder perbit --interpreter dp --cp";;
    7) echo "--method proaicd --encoder xor --interpreter dp";;
    8) echo "--method proaicd --encoder xor --interpreter dp --cp";; esac; }

# refactor one variable with one method, then retrieve every bound from scratch in
# its own process, as in the paper's measurement
run_method() {
    local m=$1 var=$2 out=$3 tmp=$4
    ${R} "${D}" "${RF}" --vars ${var} --dims ${DIMS} --dtype d $(refactor_opts $m) > "$tmp" 2>&1
    echo "Method#$m, $(method_name $m), $(grep -m1 'Refactor time' "$tmp")" >> "$out"
    for eb in ${EBS}; do
        ${Q} "${D}" "${RF}" --vars ${var} --dtype d $(retrieve_opts $m) --tolerance $eb > "$tmp" 2>&1
        echo "Method#$m, $(method_name $m), ErrorBound=$eb, $(grep -m1 'Bitrate' "$tmp"), $(grep -m1 'PSNR' "$tmp"), $(grep -m1 'Max error' "$tmp"), $(grep -m1 'Reconstruct time' "$tmp")" >> "$out"
    done
}

stage_ablation() {
    mkdir -p "${RESULT_DIR}" "${RF}"
    for var in ${VARS}; do
        out=${RESULT_DIR}/${DATASET}_${var}_ablation.txt; tmp=${RESULT_DIR}/${DATASET}_${var}_tmp.txt
        [ -f "${D}/${var}.dat" ] || { echo "missing ${D}/${var}.dat (run the prepare stage)"; exit 1; }
        > "$out"
        for m in 2 3 4 5 6 7 8; do run_method $m ${var} "$out" "$tmp"; done
        rm -f "$tmp"
        echo "${DATASET} ${var}: results in $out"
    done
}

# Figure 12: the CESM temperature field reconstructed at the paper's bounds
stage_visualize() {
    [ "${DATASET}" = "CESM" ] || { echo "the visualize stage is for DATASET=CESM"; exit 1; }
    local var=T out=${RESULT_DIR}/CESM_T_visualization.txt tmp=${RESULT_DIR}/CESM_T_vis_tmp.txt
    mkdir -p "${RESULT_DIR}" "${RF}"; > "$out"
    # method, relative bound, output name
    for spec in "3 4e-3 PMGARD_T" "4 1e-2 SZ3R_T_low" "4 1e-3 SZ3R_T_high" "8 2e-3 TWO_T"; do
        set -- $spec; m=$1; eb=$2; name=$3
        ${R} "${D}" "${RF}" --vars ${var} --dims ${DIMS} --dtype d $(refactor_opts $m) > "$tmp" 2>&1
        ${Q} "${D}" "${RF}" --vars ${var} --dtype d $(retrieve_opts $m) --tolerance $eb --output "${RF}/out_${name}" > "$tmp" 2>&1
        mv "${RF}/out_${name}/${var}.dat" "${RF}/${var}_refactored/${name}.dat" && rmdir "${RF}/out_${name}"
        echo "Method#$m, $(method_name $m), ErrorBound=$eb, $(grep -m1 'Bitrate' "$tmp"), $(grep -m1 'PSNR' "$tmp"), $(grep -m1 'Reconstruct time' "$tmp"), file=${RF}/${var}_refactored/${name}.dat" >> "$out"
    done
    rm -f "$tmp"; echo "reconstructions listed in $out"
}

# Tables 3 and 4: averages over the variables of the dataset
stage_summary() {
    python3 - "${RESULT_DIR}" "${DATASET}" ${VARS} <<'PY'
import sys, re, os
result_dir, dataset, vars_ = sys.argv[1], sys.argv[2], sys.argv[3:]
target = 1e-4
ref, rec, br = {}, {}, {}
for v in vars_:
    p = f"{result_dir}/{dataset}_{v}_ablation.txt"
    if not os.path.exists(p): print(f"missing {p}"); continue
    for l in open(p):
        m = re.match(r"Method#(\d+), ([^,]+), (.*)", l.strip())
        if not m: continue
        key = (int(m.group(1)), m.group(2)); rest = m.group(3)
        t = re.search(r"Refactor time: ([\d.]+)s", rest)
        if t: ref.setdefault(key, []).append(float(t.group(1))); continue
        eb = re.search(r"ErrorBound=([\d.eE+-]+)", rest)
        if not eb or abs(float(eb.group(1)) - target) > 0.01 * target: continue
        t = re.search(r"Reconstruct time: ([\d.]+)s", rest); b = re.search(r"Bitrate = ([\d.]+)", rest)
        if t: rec.setdefault(key, []).append(float(t.group(1)))
        if b: br.setdefault(key, []).append(float(b.group(1)))
avg = lambda d, k: (sum(d[k]) / len(d[k])) if k in d and d[k] else float("nan")
print(f"{dataset}: averages over {len(vars_)} variables")
print(f"{'method':22s}{'refactor time (s)':>20s}{'retrieval time at 1e-4 (s)':>28s}{'bitrate at 1e-4':>18s}")
for key in sorted(set(ref) | set(rec)):
    print(f"#{key[0]} {key[1]:18s}{avg(ref, key):20.3f}{avg(rec, key):28.3f}{avg(br, key):18.4f}")
PY
}

# SDRBench download and conversion to the per-variable double-precision layout
stage_prepare() {
    python3 - "${DATASET}" "${DATA_ROOT}" <<'PY'
import sys, os, glob, tarfile, urllib.request, shutil
import numpy as np
dataset, root = sys.argv[1], sys.argv[2]
BASE = "https://g-8d6b0.fd635.8443.data.globus.org/ds131.2/Data-Reduction-Repo/raw-data"
cfg = {
    "CESM":    (f"{BASE}/CESM-ATM/SDRBENCH-CESM-ATM-26x1800x3600.tar.gz", "SDRBENCH-CESM-ATM-26x1800x3600", "*.f32", np.float32, lambda f: f.split("_")[0]),
    "Miranda": (f"{BASE}/Miranda/SDRBENCH-Miranda-256x384x384.tar.gz", "SDRBENCH-Miranda-256x384x384", "*.d64", np.float64, lambda f: f.split(".")[0]),
    "SCALE":   (f"{BASE}/SCALE_LETKF/SDRBENCH-SCALE-98x1200x1200.tar.gz", "SDRBENCH-SCALE_98x1200x1200", "*.f32", np.float32, lambda f: f.split("-")[0]),
    "S3D":     (f"{BASE}/S3D/SDRBENCH-S3D.tar.gz", "SDRBENCH-S3D", None, np.float64, None),
}
if dataset not in cfg: sys.exit(f"no download recipe for {dataset}")
url, subdir, pattern, dtype, varname = cfg[dataset]
base = os.path.join(root, dataset); data_dir = os.path.join(base, "data")
os.makedirs(data_dir, exist_ok=True); os.makedirs(os.path.join(base, "refactor"), exist_ok=True)
archive = os.path.join(base, os.path.basename(url))
if not os.path.exists(archive):
    print(f"downloading {url}"); urllib.request.urlretrieve(url, archive)
print(f"extracting {archive}")
with tarfile.open(archive, "r:gz") as tar: tar.extractall(path=base)
src = os.path.join(base, subdir)
if dataset == "S3D":
    # one multi-variable file (11 x 500^3, double); variable order per SDRBENCH-S3D/template.txt
    names = ["CH4", "O2", "CO", "CO2", "H2O", "N2", "Temperature", "Pressure", "VelocityX", "VelocityY", "VelocityZ"]
    data = np.fromfile(os.path.join(src, "stat_planar.2.9950E-03.field.d64"), dtype=np.float64).reshape((11, 500, 500, 500))
    for i, n in enumerate(names):
        data[i].tofile(os.path.join(data_dir, f"{n}.dat")); print(f"  {n}.dat")
    del data
else:
    for path in sorted(glob.glob(os.path.join(src, pattern))):
        n = varname(os.path.basename(path))
        np.fromfile(path, dtype=dtype).astype(np.float64).tofile(os.path.join(data_dir, f"{n}.dat")); print(f"  {n}.dat")
shutil.rmtree(src); os.remove(archive)
print(f"{dataset} prepared under {data_dir}")
PY
}

stages=${@:-ablation}
for stage in ${stages}; do
    case $stage in
        prepare)   stage_prepare ;;
        ablation)  stage_ablation ;;
        summary)   stage_summary ;;
        visualize) stage_visualize ;;
        *) echo "unknown stage $stage (prepare|ablation|summary|visualize)"; exit 1 ;;
    esac
done
