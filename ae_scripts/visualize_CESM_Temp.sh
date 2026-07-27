#!/bin/bash
#SBATCH --time=03:00:00                 # Time limit for the job (adjust as needed)
#SBATCH --job-name=CESM_visualization   # Job name
#SBATCH --nodes=1                       # Number of nodes to allocate
#SBATCH --ntasks=16                     # Number of cores to allocate
#SBATCH --partition=                    # Partition/queue to run the job in
#SBATCH -e slurm-%j.err                 # Error file for this job
#SBATCH -o slurm-%j.out                 # Output file for this job
#SBATCH -A                              # Project allocation account name
#SBATCH --mail-type=end                 # Send email when job starts/ends
#SBATCH --mail-user=                    # Email address for notifications

data_dict_path=$1
variable="T"

# CESM 26 1800 3600

###################################################################################################
# 26 1800 3600

output_file="./ae_results/CESM_Results/CESM_${variable}_visualization.txt"
tmp_file="./ae_results/CESM_Results/CESM_${variable}_tmp.txt"

> $output_file
> $tmp_file

data_file="${data_dict_path}/data/${variable}.dat"
refactor_file="${data_dict_path}/refactor/${variable}_refactored"

###################################################################################################
# Baselines
# 1. IPComp # No need to refactor again
# ./external/IPComp/build/src/refactor $data_file "-d" "-3" 26 1800 3600 $refactor_file > $tmp_file
# time=$(grep "time" $tmp_file | head -n 1)
# echo "Method#1, IPComp, Refactor time: $time" >> $output_file

reconstructed_file="${data_dict_path}/refactor/${variable}_refactored/IPComp_T.dat"
error_bound="2e-3"

./external/IPComp/build/src/reconstructor $data_file "-d" "-3" 26 1800 3600 "-1" $error_bound $refactor_file $reconstructed_file > $tmp_file
PSNR=$(grep "PSNR" $tmp_file | head -n 1)
bitrate=$(grep "Bitrate" $tmp_file | head -n 1)
target_error=$(grep "Target" $tmp_file | head -n 1)
max_absolute_error=$(grep "Error" $tmp_file | head -n 1)
time=$(grep "time" $tmp_file | head -n 1)
echo "Method#1, IPComp, ErrorBound=$error_bound, $target_error, $max_absolute_error, $bitrate, $PSNR, $time" >> $output_file


# 3. PMGARD psnr mode
./build/test/test_mdr_refactor $data_file 3 60 3 26 1800 3600 $refactor_file 0 > $tmp_file
# time=$(grep "time" $tmp_file | head -n 1)
# echo "Method#3, PMGARD-PSNR, $time" >> $output_file

reconstructed_file="${data_dict_path}/refactor/${variable}_refactored/PMGARD_T.dat"
error_bound="4e-3"

./build/test/test_mdr_reconstructor $data_file 1 $error_bound $refactor_file 0 $reconstructed_file > $tmp_file
PSNR=$(grep "PSNR" $tmp_file | head -n 1)
bitrate=$(grep "Bitrate" $tmp_file | head -n 1)
time=$(grep "time" $tmp_file | head -n 1)
echo "Method#3, PMGARD-PSNR, ErrorBound=$error_bound, $bitrate, $PSNR, $time" >> $output_file


# 4. SZ3-R # No need to refactor again
# ./build/test/test_PSZ3-delta_refactor $data_file $refactor_file 18 3 26 1800 3600 "-d" > $tmp_file
# time=$(grep "time" $tmp_file | head -n 1)
# echo "Method#4, SZ3-R, $time" >> $output_file

reconstructed_file="${data_dict_path}/refactor/${variable}_refactored/SZ3R_T_low.dat"
error_bound="1e-2"

./build/test/test_PSZ3-delta_reconstructor $data_file $refactor_file 1 $error_bound "-d" $reconstructed_file > $tmp_file
PSNR=$(grep "PSNR" $tmp_file | head -n 1)
bitrate=$(grep "Bitrate" $tmp_file | head -n 1)
time=$(grep "time" $tmp_file | head -n 1)
echo "Method#4, SZ3-R, ErrorBound=$error_bound, $bitrate, $PSNR, $time" >> $output_file

reconstructed_file="${data_dict_path}/refactor/${variable}_refactored/SZ3R_T_high.dat"
error_bound="1e-3"

./build/test/test_PSZ3-delta_reconstructor $data_file $refactor_file 1 $error_bound "-d" $reconstructed_file > $tmp_file
PSNR=$(grep "PSNR" $tmp_file | head -n 1)
bitrate=$(grep "Bitrate" $tmp_file | head -n 1)
time=$(grep "time" $tmp_file | head -n 1)
echo "Method#4, SZ3-R, ErrorBound=$error_bound, $bitrate, $PSNR, $time" >> $output_file


# 8. CoeffDecom psnr mode
./build/test/two_modes_refactor $data_file "-d" 3 60 3 26 1800 3600 $refactor_file "-XOR" "-PSNR" "-CP" > $tmp_file
# time=$(grep "Refactor_time" $tmp_file | head -n 1)
# echo "Method#8, CoeffDecom-PSNR, $time" >> $output_file

reconstructed_file="${data_dict_path}/refactor/${variable}_refactored/TWO_T.dat"
error_bound="2e-3"

./build/test/two_modes_reconstructor $data_file "-d" 1 $error_bound $refactor_file "-XOR" "-DP" "-CP" $reconstructed_file > $tmp_file
PSNR=$(grep "PSNR" $tmp_file | head -n 1)
bitrate=$(grep "Bitrate" $tmp_file | head -n 1)
time=$(grep "Reconstruct_time" $tmp_file | head -n 1)
echo "Method#8, CoeffDecom-PSNR, ErrorBound=$error_bound, $bitrate, $PSNR, $time" >> $output_file

