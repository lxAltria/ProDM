#!/bin/bash
# ============================================================================
# GE reproduction (QPro-N / QPro-NW / QPro, 6 QoIs)
#
# Requirements:
#  - ProDM built via build_script.sh (all PRODM_WITH_* options ON; the GE QoI
#    tools qoi_{T,C,Mach,PT,mu}_d64 are built under build/app/GE/sequential
#    and refactor_d64 / qoi_Vtot_d64 under build/test)
#  - GE data under ${DATA_DIR}/data (VelocityX/Y/Z.dat, Pressure.dat,
#    Density.dat, block_sizes.dat); the refactor directory is created
#    automatically under ${DATA_DIR}/refactor
#
# Defaults assume the data was placed under build/app/GE/data; override via
# environment, e.g.:
#   DATA_DIR=/path/to/GE_small_reordered sh GE_steps.sh
# ============================================================================

script_dir=$(cd "$(dirname "$0")" && pwd)
PRODM_DIR=${PRODM_DIR:-$(cd "${script_dir}/../.." && pwd)}
build_dir=${PRODM_DIR}/build
DATA_DIR=${DATA_DIR:-${build_dir}/app/GE}
RESULT_DIR=${RESULT_DIR:-${build_dir}/Result}

refactor_cmd=${build_dir}/test/refactor_d64

# qoi_Vtot_d64 is shared with the example workflow and lives under test/;
# the GE-specific QoI tools live under app/GE/sequential/
qoi_cmd() {
    if [ "$1" = "Vtot" ]; then
        echo "${build_dir}/test/qoi_Vtot_d64"
    else
        echo "${build_dir}/app/GE/sequential/qoi_$1_d64"
    fi
}

data="GE"
QoIs=("Vtot" "T" "C" "Mach" "PT" "mu")

# error bounds: two geometric ladders (0.1 and 0.05, ratio 0.1, 6 steps each)
a1=0.1
a2=0.05
r=0.1
n=6
error_bounds=()
a=$a1
for ((i = 1; i <= n; i++)); do
    error_bounds+=($a)
    a=$(echo "scale=10; $a * $r" | bc)
done
a=$a2
for ((i = 1; i <= n; i++)); do
    error_bounds+=($a)
    a=$(echo "scale=10; $a * $r" | bc)
done
error_bounds=($(printf "%s\n" "${error_bounds[@]}" | sort -nr))

# GEApproximator expects block_sizes.dat at the dataset root (next to data/ and refactor/)
if [ ! -f "${DATA_DIR}/block_sizes.dat" ] && [ -f "${DATA_DIR}/data/block_sizes.dat" ]; then
    cp "${DATA_DIR}/data/block_sizes.dat" "${DATA_DIR}/block_sizes.dat"
fi

mkdir -p ${RESULT_DIR}
output_file=${RESULT_DIR}/GE_result_3D.txt
tmp_file=${RESULT_DIR}/GE_tmp.txt
> $output_file
> $tmp_file

run_retrievals() {
    method=$1; weighted=$2; decrease_method=$3
    for qoi in "${QoIs[@]}"; do
        cmd=$(qoi_cmd $qoi)
        for error_bound in "${error_bounds[@]}"; do
            $cmd 3 $weighted $decrease_method "$error_bound" "${DATA_DIR}" > $tmp_file
            bitrate=$(grep "bitrate" $tmp_file | head -n 1)
            time=$(grep "elapsed_time" $tmp_file | head -n 1)
            requested_tau=$(grep "requested_error" $tmp_file | head -n 1)
            max_est_err=$(grep "max_est_error" $tmp_file | head -n 1)
            max_act_err=$(grep "max_act_error" $tmp_file | head -n 1)
            echo "$data, Method#$method: $qoi, ErrorBound=$error_bound, $time, $requested_tau, $max_est_err, $max_act_err, $bitrate" >> $output_file
        done
    done
}

# 4. QPro-N (unweighted)
${refactor_cmd} 3 0 $data "${DATA_DIR}" > $tmp_file
time=$(grep "elapsed_time" $tmp_file | head -n 1)
echo "$data, Method#4: Refactor: $time" >> $output_file
run_retrievals 4 0 0

# 5. QPro-NW (weighted refactor, naive eb decrease)
${refactor_cmd} 3 1 $data "${DATA_DIR}" 4 3 0.001 > $tmp_file
time=$(grep "elapsed_time" $tmp_file | head -n 1)
echo "$data, Method#5: Refactor: $time" >> $output_file
run_retrievals 5 1 0

# 6. QPro (weighted refactor, coordinate eb decrease)
${refactor_cmd} 3 1 $data "${DATA_DIR}" 4 3 0.001 > $tmp_file
time=$(grep "elapsed_time" $tmp_file | head -n 1)
echo "$data, Method#6: Refactor: $time" >> $output_file
run_retrievals 6 1 1

rm $tmp_file
echo "GE reproduction finished; results in $output_file"
