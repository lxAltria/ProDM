#!/bin/bash
# ============================================================================
# HPDC'26 (QProR): QProR-N / QProR-NW / QProR on a structured dataset (Hurricane, NYX, SCALE, S3D), QoI V_total.
#
# Requirements: ProDM built via build_script.sh (the GE tools refactor_d64 and qoi_Vtot_d64 under build/app/GE);
# the dataset arranged as ${DATA_DIR}/data/{VelocityX,VelocityY,VelocityZ}.dat (double precision; the grid dimensions
# are known to the tools from the dataset name). The refactored representation is written to ${DATA_DIR}/refactor.
#
# usage: DATASET=Hurricane DATA_DIR=/path/to/Hurricane_d64 sh structured_steps.sh
# optional: APPROX (approximator id, default 1 = SZ3-based approximation; 4 = HPEZ), MAX_WEIGHT (default 7),
#           BLOCK_SIZE (default 4), APPROX_EB (default 0.001), RESULT_DIR
# ============================================================================
set -u
script_dir=$(cd "$(dirname "$0")" && pwd)
PRODM_DIR=${PRODM_DIR:-$(cd "${script_dir}/../.." && pwd)}
BIN=${PRODM_DIR}/build/app/GE
DATASET=${DATASET:?set DATASET to Hurricane, NYX, SCALE or S3D}
DATA_DIR=${DATA_DIR:?set DATA_DIR to the dataset directory (holding data/ )}
APPROX=${APPROX:-1}; MAX_WEIGHT=${MAX_WEIGHT:-7}; BLOCK_SIZE=${BLOCK_SIZE:-4}; APPROX_EB=${APPROX_EB:-0.001}
RESULT_DIR=${RESULT_DIR:-${PRODM_DIR}/build/Result}
error_bounds=${EBS:-"0.1 0.05 0.01 0.005 0.001 0.0005 0.0001 0.00005 0.00001 0.000005 0.000001 0.0000005"}
mkdir -p ${RESULT_DIR}
out=${RESULT_DIR}/${DATASET}_result.txt; tmp=${RESULT_DIR}/${DATASET}_tmp.txt
> $out
run_retrievals() {  # method_name weighted decrease_method
    for eb in $error_bounds; do
        ${BIN}/qoi_Vtot_d64 ${APPROX} $2 $3 $eb ${DATA_DIR} > $tmp
        echo "${DATASET}, $1: Vtot, ErrorBound=$eb, $(grep -m1 elapsed_time $tmp), $(grep -m1 requested_error $tmp), $(grep -m1 max_est_error $tmp), $(grep -m1 max_act_error $tmp), $(grep -m1 bitrate $tmp)" >> $out
    done
}
# QProR-N: plain bitplane encoding, uniform tightening
${BIN}/refactor_d64 ${APPROX} 0 ${DATASET} ${DATA_DIR} > $tmp
echo "${DATASET}, QProR-N: Refactor: $(grep -m1 elapsed_time $tmp)" >> $out
run_retrievals "QProR-N" 0 0
# QProR-NW: weighted bitplane encoding, uniform tightening; QProR: weighted encoding, coordinate error-bound descent
${BIN}/refactor_d64 ${APPROX} 1 ${DATASET} ${DATA_DIR} ${MAX_WEIGHT} ${BLOCK_SIZE} ${APPROX_EB} > $tmp
echo "${DATASET}, QProR-NW/QProR: Refactor: $(grep -m1 elapsed_time $tmp)" >> $out
run_retrievals "QProR-NW" 1 0
run_retrievals "QProR" 1 1
rm -f $tmp
echo "done; results in $out"
