#!/bin/bash
#SBATCH --time=03:00:00               # Time limit for the job (adjust as needed)
#SBATCH --job-name=weak_scale_JHTDB  # Job name
#SBATCH --ntasks=1024
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12G
#SBATCH --partition=                  # Partition/queue to run the job in
#SBATCH -e slurm-%j.err               # Error file for this job
#SBATCH -o slurm-%j.out               # Output file for this job
#SBATCH -A                            # Project allocation account name
#SBATCH --mail-type=end               # Send email when job starts/ends
#SBATCH --mail-user=                  # Email address for notifications

data_dict_path=$1

output_file="./ae_results/JHTDB_Results/results_weak_scale.txt"
tmp_file="./ae_results/JHTDB_Results/results_tmp.txt"

> $output_file
> $tmp_file

data_file="${data_dict_path}/data/Temperature_"
refactor_file="${data_dict_path}/"

error_bound=1e-4
core_counts=(128 256 512 1024)
NUM_RUNS=1

# Helper function: extract elapsed_time value from a line like "elapsed_time: 1.234"
extract_time() {
    grep "elapsed_time" "$1" | head -n 1 | grep -oE '[0-9]+\.?[0-9]*'
}

# Helper function: compute average using awk
compute_avg() {
    local -n arr=$1
    printf '%s\n' "${arr[@]}" | awk '{s+=$1} END{printf "%.6f", s/NR}'
}

for np in "${core_counts[@]}"; do
    echo "========== Cores: $np, ErrorBound: $error_bound ==========" >> $output_file

    # ======== Method#1: PMGARD-HB ========
    write_file="${data_dict_path}/PMGARD/$np"

    refactor_times=()
    recon_times=()
    for run in $(seq 1 $NUM_RUNS); do
        mpirun -n $np ./build/prl_src/para_mdr_refactor $data_file 4 60 3 256 512 512 $refactor_file > $tmp_file
        t=$(extract_time $tmp_file)
        refactor_times+=("$t")

        mpirun -n $np ./build/prl_src/para_mdr_reconstructor $data_file 1 $error_bound $refactor_file $write_file > $tmp_file
        t=$(extract_time $tmp_file)
        recon_times+=("$t")

        # Save last run's metrics
        if [ $run -eq $NUM_RUNS ]; then
            bitrate=$(grep "bitrate" $tmp_file | head -n 1)
            requested_tau=$(grep "tolerance" $tmp_file | head -n 1)
            max_act_err=$(grep "max_abs_error" $tmp_file | head -n 1)
        fi
    done
    avg_refactor=$(compute_avg refactor_times)
    avg_recon=$(compute_avg recon_times)
    echo "Method#1, Cores=$np, Refactor: avg_elapsed_time=$avg_refactor (runs: ${refactor_times[*]})" >> $output_file
    echo "Method#1, Cores=$np, ErrorBound=$error_bound, avg_elapsed_time=$avg_recon (runs: ${recon_times[*]}), $requested_tau, $max_act_err, $bitrate" >> $output_file

    # ======== Method#2: SZ3-R ========
    write_file="${data_dict_path}/SZ3-R/$np"

    refactor_times=()
    recon_times=()
    for run in $(seq 1 $NUM_RUNS); do
        mpirun -n $np ./build/prl_src/para_PSZ3-delta_refactor $data_file $refactor_file 18 3 256 512 512 "-d" > $tmp_file
        t=$(extract_time $tmp_file)
        refactor_times+=("$t")

        mpirun -n $np ./build/prl_src/para_PSZ3-delta_reconstructor $data_file $refactor_file 1 $error_bound "-d" $write_file > $tmp_file
        t=$(extract_time $tmp_file)
        recon_times+=("$t")

        if [ $run -eq $NUM_RUNS ]; then
            bitrate=$(grep "bitrate" $tmp_file | head -n 1)
            requested_tau=$(grep "tolerance" $tmp_file | head -n 1)
            max_act_err=$(grep "max_abs_error" $tmp_file | head -n 1)
        fi
    done
    avg_refactor=$(compute_avg refactor_times)
    avg_recon=$(compute_avg recon_times)
    echo "Method#2, Cores=$np, Refactor: avg_elapsed_time=$avg_refactor (runs: ${refactor_times[*]})" >> $output_file
    echo "Method#2, Cores=$np, ErrorBound=$error_bound, avg_elapsed_time=$avg_recon (runs: ${recon_times[*]}), $requested_tau, $max_act_err, $bitrate" >> $output_file

    # ======== Method#3: IPComp ========
    write_file="${data_dict_path}/IPComp/$np"

    refactor_times=()
    recon_times=()
    for run in $(seq 1 $NUM_RUNS); do
        mpirun -n $np ./external/IPComp/build/prl_src/para_IPComp_refactor $data_file "-d" "-3" 256 512 512 $refactor_file > $tmp_file
        t=$(extract_time $tmp_file)
        refactor_times+=("$t")

        mpirun -n $np ./external/IPComp/build/prl_src/para_IPComp_reconstructor $data_file "-d" "-3" 256 512 512 "-1" $error_bound $refactor_file $write_file > $tmp_file
        t=$(extract_time $tmp_file)
        recon_times+=("$t")

        if [ $run -eq $NUM_RUNS ]; then
            bitrate=$(grep "bitrate" $tmp_file | head -n 1)
            requested_tau=$(grep "tolerance" $tmp_file | head -n 1)
            max_act_err=$(grep "max_abs_error" $tmp_file | head -n 1)
        fi
    done
    avg_refactor=$(compute_avg refactor_times)
    avg_recon=$(compute_avg recon_times)
    echo "Method#3, Cores=$np, Refactor: avg_elapsed_time=$avg_refactor (runs: ${refactor_times[*]})" >> $output_file
    echo "Method#3, Cores=$np, ErrorBound=$error_bound, avg_elapsed_time=$avg_recon (runs: ${recon_times[*]}), $requested_tau, $max_act_err, $bitrate" >> $output_file

    # ======== Method#4: TWO - eb ========
    write_file="${data_dict_path}/TWO/$np"

    refactor_times=()
    recon_times=()
    for run in $(seq 1 $NUM_RUNS); do
        mpirun -n $np ./build/prl_src/para_two_modes_refactor $data_file "-d" 4 60 3 256 512 512 $refactor_file "-PerBit" "-eb" "-CP" > $tmp_file
        t=$(extract_time $tmp_file)
        refactor_times+=("$t")

        mpirun -n $np ./build/prl_src/para_two_modes_reconstructor $data_file "-d" 1 $error_bound $refactor_file "-PerBit" "-DP" "-CP" $write_file > $tmp_file
        t=$(extract_time $tmp_file)
        recon_times+=("$t")

        if [ $run -eq $NUM_RUNS ]; then
            bitrate=$(grep "bitrate" $tmp_file | head -n 1)
            requested_tau=$(grep "tolerance" $tmp_file | head -n 1)
            max_act_err=$(grep "max_abs_error" $tmp_file | head -n 1)
        fi
    done
    avg_refactor=$(compute_avg refactor_times)
    avg_recon=$(compute_avg recon_times)
    echo "Method#4, Cores=$np, Refactor: avg_elapsed_time=$avg_refactor (runs: ${refactor_times[*]})" >> $output_file
    echo "Method#4, Cores=$np, ErrorBound=$error_bound, avg_elapsed_time=$avg_recon (runs: ${recon_times[*]}), $requested_tau, $max_act_err, $bitrate" >> $output_file
done

rm $tmp_file