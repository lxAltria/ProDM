#!/bin/bash
#SBATCH --time=03:00:00                 # Time limit for the job (adjust as needed)
#SBATCH --job-name=CESM_ablation_study   # Job name
#SBATCH --nodes=1                       # Number of nodes to allocate
#SBATCH --ntasks=16                     # Number of cores to allocate
#SBATCH --partition=                    # Partition/queue to run the job in
#SBATCH -e slurm-%j.err                 # Error file for this job
#SBATCH -o slurm-%j.out                 # Output file for this job
#SBATCH -A                              # Project allocation account name
#SBATCH --mail-type=end                 # Send email when job starts/ends
#SBATCH --mail-user=                    # Email address for notifications

error_bounds=(0.1 0.05 0.01 0.005 0.001 0.0005 0.0001 0.00005 0.00001 0.000005 0.000001 0.0000005 0.0000001 0.00000005 0.00000001 0.000000005 0.000000001)

data_dict_path=$1
variable=$2

# CESM 26 1800 3600

###################################################################################################
# 26 1800 3600
output_file="./ae_results/CESM_Results/CESM_${variable}_ablation.txt"
tmp_file="./ae_results/CESM_Results/CESM_${variable}_tmp.txt"

> $output_file
> $tmp_file

data_file="${data_dict_path}/data/${variable}.dat"
refactor_file="${data_dict_path}/refactor/${variable}_refactored"

###################################################################################################
# Baselines
# 1. IPComp
./external/IPComp/build/src/refactor $data_file "-d" "-3" 26 1800 3600 $refactor_file > $tmp_file
time=$(grep "time" $tmp_file | head -n 1)
echo "Method#1, IPComp, Refactor time: $time" >> $output_file
for error_bound in "${error_bounds[@]}"; do
    ./external/IPComp/build/src/reconstructor $data_file "-d" "-3" 26 1800 3600 "-1" $error_bound $refactor_file > $tmp_file
    PSNR=$(grep "PSNR" $tmp_file | head -n 1)
    bitrate=$(grep "Bitrate" $tmp_file | head -n 1)
    target_error=$(grep "Target" $tmp_file | head -n 1)
    max_absolute_error=$(grep "Error" $tmp_file | head -n 1)
    time=$(grep "time" $tmp_file | head -n 1)
    echo "Method#1, IPComp, ErrorBound=$error_bound, $target_error, $max_absolute_error, $bitrate, $PSNR, $time" >> $output_file
done

# 2. PMGARD eb mode
./build/test/test_mdr_refactor $data_file 3 60 3 26 1800 3600 $refactor_file 1 > $tmp_file
time=$(grep "time" $tmp_file | head -n 1)
echo "Method#2, PMGARD-EB, $time" >> $output_file
for error_bound in "${error_bounds[@]}"; do
    ./build/test/test_mdr_reconstructor $data_file 1 $error_bound $refactor_file 1 > $tmp_file
    PSNR=$(grep "PSNR" $tmp_file | head -n 1)
    bitrate=$(grep "Bitrate" $tmp_file | head -n 1)
    time=$(grep "time" $tmp_file | head -n 1)
    echo "Method#2, PMGARD-EB, ErrorBound=$error_bound, $bitrate, $PSNR, $time" >> $output_file
done

# 3. PMGARD psnr mode
./build/test/test_mdr_refactor $data_file 3 60 3 26 1800 3600 $refactor_file 0 > $tmp_file
time=$(grep "time" $tmp_file | head -n 1)
echo "Method#3, PMGARD-PSNR, $time" >> $output_file
for error_bound in "${error_bounds[@]}"; do
    ./build/test/test_mdr_reconstructor $data_file 1 $error_bound $refactor_file 0 > $tmp_file
    PSNR=$(grep "PSNR" $tmp_file | head -n 1)
    bitrate=$(grep "Bitrate" $tmp_file | head -n 1)
    time=$(grep "time" $tmp_file | head -n 1)
    echo "Method#3, PMGARD-PSNR, ErrorBound=$error_bound, $bitrate, $PSNR, $time" >> $output_file
done

# 4. SZ3-R
./build/test/test_PSZ3-delta_refactor $data_file $refactor_file 18 3 26 1800 3600 "-d" > $tmp_file
time=$(grep "time" $tmp_file | head -n 1)
echo "Method#4, SZ3-R, $time" >> $output_file
for error_bound in "${error_bounds[@]}"; do
    ./build/test/test_PSZ3-delta_reconstructor $data_file $refactor_file 1 $error_bound "-d" > $tmp_file
    PSNR=$(grep "PSNR" $tmp_file | head -n 1)
    bitrate=$(grep "Bitrate" $tmp_file | head -n 1)
    time=$(grep "time" $tmp_file | head -n 1)
    echo "Method#4, SZ3-R, ErrorBound=$error_bound, $bitrate, $PSNR, $time" >> $output_file
done

###################################################################################################
# Prior-eb
# 5. AdatInterp eb mode
./build/test/two_modes_refactor $data_file "-d" 3 60 3 26 1800 3600 $refactor_file "-PerBit" "-eb" "-no_CP" > $tmp_file
time=$(grep "Refactor_time" $tmp_file | head -n 1)
echo "Method#5, AdatInterp-EB, $time" >> $output_file
for error_bound in "${error_bounds[@]}"; do
    ./build/test/two_modes_reconstructor $data_file "-d" 1 $error_bound $refactor_file "-PerBit" "-DP" "-no_CP"  > $tmp_file
    PSNR=$(grep "PSNR" $tmp_file | head -n 1)
    bitrate=$(grep "Bitrate" $tmp_file | head -n 1)
    time=$(grep "Reconstruct_time" $tmp_file | head -n 1)
    echo "Method#5, AdatInterp-EB, ErrorBound=$error_bound, $bitrate, $PSNR, $time" >> $output_file
done

# 6. CoeffDecom eb mode
./build/test/two_modes_refactor $data_file "-d" 3 60 3 26 1800 3600 $refactor_file "-PerBit" "-eb" "-CP" > $tmp_file
time=$(grep "Refactor_time" $tmp_file | head -n 1)
echo "Method#6, CoeffDecom-EB, $time" >> $output_file
for error_bound in "${error_bounds[@]}"; do
    ./build/test/two_modes_reconstructor $data_file "-d" 1 $error_bound $refactor_file "-PerBit" "-DP" "-CP" > $tmp_file
    PSNR=$(grep "PSNR" $tmp_file | head -n 1)
    bitrate=$(grep "Bitrate" $tmp_file | head -n 1)
    time=$(grep "Reconstruct_time" $tmp_file | head -n 1)
    echo "Method#6, CoeffDecom-EB, ErrorBound=$error_bound, $bitrate, $PSNR, $time" >> $output_file
done

###################################################################################################
# Prior-PSNR
# 7. AdatInterp psnr mode
./build/test/two_modes_refactor $data_file "-d" 3 60 3 26 1800 3600 $refactor_file "-XOR" "-PSNR" "-no_CP" > $tmp_file
time=$(grep "Refactor_time" $tmp_file | head -n 1)
echo "Method#7, AdatInterp-PSNR, $time" >> $output_file
for error_bound in "${error_bounds[@]}"; do
    ./build/test/two_modes_reconstructor $data_file "-d" 1 $error_bound $refactor_file "-XOR" "-DP" "-no_CP"  > $tmp_file
    PSNR=$(grep "PSNR" $tmp_file | head -n 1)
    bitrate=$(grep "Bitrate" $tmp_file | head -n 1)
    time=$(grep "Reconstruct_time" $tmp_file | head -n 1)
    echo "Method#7, AdatInterp-PSNR, ErrorBound=$error_bound, $bitrate, $PSNR, $time" >> $output_file
done

# 8. CoeffDecom psnr mode
./build/test/two_modes_refactor $data_file "-d" 3 60 3 26 1800 3600 $refactor_file "-XOR" "-PSNR" "-CP" > $tmp_file
time=$(grep "Refactor_time" $tmp_file | head -n 1)
echo "Method#8, CoeffDecom-PSNR, $time" >> $output_file
for error_bound in "${error_bounds[@]}"; do
    ./build/test/two_modes_reconstructor $data_file "-d" 1 $error_bound $refactor_file "-XOR" "-DP" "-CP" > $tmp_file
    PSNR=$(grep "PSNR" $tmp_file | head -n 1)
    bitrate=$(grep "Bitrate" $tmp_file | head -n 1)
    time=$(grep "Reconstruct_time" $tmp_file | head -n 1)
    echo "Method#8, CoeffDecom-PSNR, ErrorBound=$error_bound, $bitrate, $PSNR, $time" >> $output_file
done

rm $tmp_file
